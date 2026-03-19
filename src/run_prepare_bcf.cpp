#include <htslib/vcf.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include "prepare_bcf.h"
#include <boost/program_options.hpp>

using namespace std;
namespace po=boost::program_options;

void mergeWithHetFile ( 
    const bcf_hdr_t *hdr, 
    const bcf1_t *rec, 
    variantAD_t *lastVar, 
    ifstream *hetFiles, 
    int *GT, 
    int *AD,
    bool genoQC,
    int *numDiscarded
) {
    for (int i = 0; i < bcf_hdr_nsamples(hdr); i++) {

        if (lastVar[i].rid < rec->rid) {
            // the variant in the HET file is on an earlier chom
            lastVar[i] = consumeVariants(&hetFiles[i], rec);
        } else if (lastVar[i].pos < rec->pos && lastVar[i].rid == rec->rid) {
            // the variant in the HET file comes before the current imputed variant
            lastVar[i] = consumeVariants(&hetFiles[i], rec);
        }
        if (checkVariantEquality(lastVar[i], rec) != 1) {
            // if the variant in the HET file comes after the current imputed variant
            if (bcf_gt_allele(GT[2*i]) + bcf_gt_allele(GT[2*i+1]) == 1 && genoQC) {
                // set genotype to missing if imputed as a het
                GT[2*i] = GT[2*i+1] = bcf_gt_missing;
                numDiscarded[i]+=(rec->rid+1 <= 22);
            }
            continue;
        }
        int discarded = writeADs(i, AD, lastVar[i], GT);
        if (discarded && genoQC) {
            // set genotype as missing if not imputed as a het
            GT[2*i] = GT[2*i+1] = bcf_gt_missing;
            numDiscarded[i]+=(rec->rid+1 <= 22);
        }
    }
}

int prepare_bcf_main(std::vector<std::string> &argv) {

    const char *cnv_mask_fin = NULL;
    const char *ref_bias_fin = NULL;
    const char *ref_bias_fout = NULL;
    const char *hetFileFormatStr = NULL;
    bool clearADs = false, outputADs = false, genoQC = false, mask_deletions = false;

    po::variables_map vm;
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "produce help message")
        ("clear-ADs", "flag for clearing ADs already in the file")
        ("write-ADs", "flag to set whether or not ADs are written in output")
        ("geno-QC", "flag to set whether or not to QC genotypes")
        ("mask-deletions", "flag to set whether or not to mask deletions called from depth")
        ("cnv-mask", po::value<std::string>(), "CNV mask file")
        ("ref-bias-in", po::value<std::string>(), "file to read ref bias from")
        ("ref-bias-out", po::value<std::string>(), "file to write computed ref bias to")
        ("het-format-str", po::value<std::string>(), "het file format str, e.g.: gvcf/%s.c21.hets.txt\n");

	try {
		po::store(po::command_line_parser(argv).options(desc).run(), vm);
		po::notify(vm);

    } catch(exception& e) {
        cerr << "error: " << e.what() << endl;
        return 1;
    }
    catch(...) {
        cerr << "Exception of unknown type!" << endl;
        return 1;
    }

    if (vm.count("help")) {
        cerr << desc << endl;
        return 0;
    }
    if (vm.count("clear-ADs")) {
        clearADs = true;
    }
    if (vm.count("geno-QC")) {
        genoQC = true;
        cerr << "Performing QC on the genotypes" << endl;
    }
    if (vm.count("write-ADs")) {
        outputADs = true;
        cerr << "Writing updated ADs to output" << endl;
    } else {
        cerr << "ADs remain unchanged in the output" << endl;
    }
    if (vm.count("ref-bias-in")) {
        cerr << "Reading ref bias from " << vm["ref-bias-in"].as<std::string>() << endl;
        ref_bias_fin = vm["ref-bias-in"].as<std::string>().c_str();
    } else {
        cerr << "No ref bias supplied" << endl;
    }	
    if (vm.count("ref-bias-out")) {
        cerr << "Writing ref bias to " << vm["ref-bias-out"].as<std::string>() << endl;
        ref_bias_fout = vm["ref-bias-out"].as<std::string>().c_str();
    }
    if (vm.count("cnv-mask")) {
        cerr << "Masking CNVs from " << vm["cnv-mask"].as<std::string>() << endl;
        cnv_mask_fin = vm["cnv-mask"].as<std::string>().c_str();
        if (vm.count("mask-deletions")) {
            mask_deletions = true;
            cerr << "Masking duplications and deletions" << endl;
        } else {
            cerr << "Only masking called duplications" << endl;
        }
    } else {
        cerr << "No CNVs being masked" << endl;
    }
    if (vm.count("het-format-str")) {
        cerr << "Looking for het ADs in  " << vm["het-format-str"].as<std::string>() << endl;
        hetFileFormatStr = vm["het-format-str"].as<std::string>().c_str();
    } 

    htsFile *fin = hts_open("-", "r"); assert(fin != NULL);
    htsFile *fout = hts_open("-", "wb"); assert(fout != NULL);
    hts_set_threads(fout, 1); // create additional thread to help compress output bcf

    bcf_hdr_t *hdr = bcf_hdr_read(fin); assert(hdr != NULL);
    const char *hdr_ad = "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths\">";
    if (!bcf_hdr_idinfo_exists(hdr, BCF_HL_FMT, bcf_hdr_id2int(hdr, BCF_DT_ID, "AD"))) {
        if(outputADs) assert(bcf_hdr_append(hdr, hdr_ad) == 0);
        clearADs = true;
        cerr << "Warning: no AD field in the input bcf" << endl;
    }
    assert(bcf_hdr_write(fout, hdr) == 0);
    int n = bcf_hdr_nsamples(hdr);

    // populate ref bias mask
    map<string, double> refBiasMap;
    if (ref_bias_fin) {
        ifstream refBiasIn;
        refBiasIn.open(ref_bias_fin); assert(refBiasIn);
        refBiasMap = parseRefBias(&refBiasIn);
        refBiasIn.close();
    }

    ofstream refBiasOut;
    if (ref_bias_fout) {
        refBiasOut.open(ref_bias_fout); assert(refBiasOut);
    }

    // populate cnv lookup table
    map<string, queue<CNV>> cnvLookup;
    std::vector<std::queue<CNV>> cnvLookupArray(n);
    if (cnv_mask_fin) {
        ifstream germlineCNVsFile;
        germlineCNVsFile.open(cnv_mask_fin); assert(germlineCNVsFile);
        cnvLookup = parseCNVs(&germlineCNVsFile, mask_deletions);
        germlineCNVsFile.close();

        for (int i = 0; i < n; i++) {
            cerr << "Number of CNVs masked in ";
            cerr << hdr->samples[i] << ": ";
            cerr << cnvLookup[hdr->samples[i]].size() << endl;
            cnvLookupArray[i] = cnvLookup[hdr->samples[i]];
        }
        cerr << "Number of regions blacklisted: " << cnvLookup["BLACKLIST"].size() << endl;
    }

    // create file handles for all gvcf files
    int numDiscarded[n];
    memset(numDiscarded, 0, sizeof(numDiscarded));
    ifstream hetFiles[n];
    variantAD_t lastVar[n];
    if (hetFileFormatStr) {
        char filename[1000];
        for (int i = 0; i < n; i++) {
            snprintf(filename, sizeof(filename), hetFileFormatStr, hdr->samples[i]);
            hetFiles[i].open(filename); 
            if (! hetFiles[i] ) {
                cerr << "Unable to open file: " << filename << endl;
            }
            lastVar[i].pos = 0;
            lastVar[i].rid = -1;
        }
    }

    bcf1_t *rec = bcf_init();
    int *GT = NULL, GT_buf_size = 0;
    int *AD = NULL, AD_buf_size = 0;

    char variant_name[1000];

    // stream the variants in the phased bcf file
    while (bcf_read(fin, hdr, rec) == 0) {
        // blacklist regions
        if (cnv_mask_fin && variantInCNV(rec->rid+1, rec->pos+1, cnvLookup["BLACKLIST"])) {
            continue;
        }
        assert(rec->n_allele == 2);
        int nGT = bcf_get_genotypes(hdr, rec, &GT, &GT_buf_size);
        assert(nGT == 2*n);
        assert(bcf_unpack(rec, BCF_UN_STR) == 0);
        snprintf(variant_name, sizeof(variant_name), "%d_%ld_%s_%s", rec->rid+1, rec->pos+1, rec->d.allele[0], rec->d.allele[1]);

        if (!clearADs) {
            int nAD = bcf_get_format_int32(hdr, rec, "AD", &AD, &AD_buf_size);
            if (nAD != 2*n) {
                cerr << "error: failed to unpack ADs from input bcf file" << endl;
                exit(1);
            }
        } else {
            if (AD_buf_size==0) {
                AD_buf_size=2*n;
                AD = (int *) malloc(sizeof(int) * AD_buf_size);
            }
            for (int i = 0; i < n; i++) {
                AD[2*i]=AD[2*i+1]=bcf_int32_missing;
            }
        }

        if (hetFileFormatStr) mergeWithHetFile(hdr, rec, lastVar, hetFiles, GT, AD, genoQC, numDiscarded);
        bcf_update_genotypes(hdr, rec, GT, GT_buf_size);

        if (cnv_mask_fin) maskCNV(rec->rid+1, rec->pos+1, hdr, cnvLookupArray, AD);
        if (ref_bias_fin) adjustRefBias(n, (refBiasMap.count(variant_name)) ? refBiasMap[variant_name] : -1, AD);
        if (ref_bias_fout) computeRefBias(AD, n, rec, &refBiasOut);
        if (outputADs) bcf_update_format_int32(hdr, rec, "AD", AD, AD_buf_size);
        assert(bcf_write(fout, hdr, rec) == 0);
    }
    free(GT);
    free(AD);

    if (hetFileFormatStr) {
        for (int i = 0; i < n; i++) {
            if (genoQC) {
                cerr << "Number of autosome sites QC'd for " << hdr->samples[i] << ": " << numDiscarded[i] << endl;
            }
            hetFiles[i].close();
        }
    }
    if (ref_bias_fout) refBiasOut.close();

    bcf_destroy(rec);
    bcf_hdr_destroy(hdr);
    hts_close(fin);
    hts_close(fout);

    return 0;
}

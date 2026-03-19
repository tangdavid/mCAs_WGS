#include <htslib/vcf.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/program_options.hpp>

using namespace std;
namespace po=boost::program_options;
namespace bio=boost::iostreams;

map<string, float> parsePRSCoeff(const char *prs_file_name) {
    map<string, float> res;
    ifstream inFile;
    inFile.open(prs_file_name, std::ios_base::in | std::ios_base::binary); assert(inFile);
    bio::filtering_istream betasFile;
    betasFile.push(bio::gzip_decompressor());
    betasFile.push(inFile);

    string fileLine; getline(betasFile, fileLine);
    string rsID, a1, a0;
    int chrom, pos;
    float beta, genpos;
    while (betasFile >> rsID >> chrom >> pos >> genpos >> a1 >> a0 >> beta) {
        res[rsID] = -beta;
    }
    return res;
}

int compute_polygenic_drive_main(std::vector<std::string> &argv) {

    const char *prs_fin = NULL;

    po::variables_map vm;
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "produce help message")
        ("prs", po::value<std::string>(), "file with PRS coefficients\n");

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
    if (vm.count("prs")) {
        cerr << "Reading PRS coefficients from " << vm["prs"].as<std::string>() << endl;
        prs_fin = vm["prs"].as<std::string>().c_str();
    } else {
        cerr << "Error: need to supply a file with PRS coefficients" << endl;
        return 1;
    }	
    htsFile *fin = hts_open("-", "r"); assert(fin != NULL);
    bcf_hdr_t *hdr = bcf_hdr_read(fin); assert(hdr != NULL);
    int n = bcf_hdr_nsamples(hdr);

    bcf1_t *rec = bcf_init();

    map<string, float> betas = parsePRSCoeff(prs_fin);
    int *AS = NULL, AS_buf_size=0;

    float *differential_PRS = (float *) calloc(n, sizeof(float));

    int used_variants = 0, discarded_variants = 0;

    while (bcf_read(fin, hdr, rec) == 0) {
        assert(rec->n_allele == 2);
        int nAS = bcf_get_format_int32(hdr, rec, "AS", &AS, &AS_buf_size);
        // cerr << rec->pos+1 << "\t" << nAS << "\t" << n << endl;
        if(nAS != n) {
            discarded_variants += 1;
            continue;
        }
        // assert(nAS == n);

        assert(bcf_unpack(rec, BCF_UN_STR) == 0);
        string rsID = rec->d.id;
        float beta;
        if(betas.count(rsID)) {
            beta = betas[rsID];
            used_variants += 1;
        } else {
            discarded_variants += 1;
            continue;
        }

        for (int i=0; i < n; i++) {
            if(AS[i] != 1 && AS[i] != -1) { 
                continue;    
            }
            differential_PRS[i] += -AS[i] * beta;
        }
    }
    
    cout << "ID" << "\t" << "PRS_differential" << endl;
    for (int i=0; i < n; i++) {
        cout << hdr->samples[i] << "\t" << differential_PRS[i] << endl;
    }
    cerr << "SNP array variants: " << betas.size() << endl;
    cerr << "Variants used: " << used_variants << endl;
    cerr << "Variants discarded: " << discarded_variants << endl;

    free(AS);
    free(differential_PRS);

    bcf_destroy(rec);
    bcf_hdr_destroy(hdr);
    hts_close(fin);

    return 0;
}


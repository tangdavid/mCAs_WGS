#include <htslib/vcf.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <boost/program_options.hpp>
#include <ctime>
#include <set>
#include <Eigen/Dense>
#include "prepare_bcf.h"
#include "gwas_tools.h"

using namespace std;
namespace po=boost::program_options;

// mCA_phase[i] == 1 if AS[i] == -1 and GT[2*i] == 0 or AS[i] == 1 and GT[2*i] == 1; else 0
// AS[i] == -1 means that mCA_phase == 0 if GT[2*i] == 1 and mCA_phase == 1 if GT[2*i] == 0; else 1
#define MCA_PHASE_IS_1(as,first_allele) ((as == -1 && bcf_gt_allele(first_allele) == 0) || (as == 1 && bcf_gt_allele(first_allele) == 1))
#define ALLELIC_SHIFT_IS_1(phase,first_allele) ((phase == 0 && bcf_gt_allele(first_allele) == 0) || (phase == 1 && bcf_gt_allele(first_allele) == 1))

Eigen::MatrixXf compute_GRM(
    const Eigen::MatrixXf &genos
) {
    // Compute the GRM using matrix multiplication
    float m = genos.rows();
    
    Eigen::MatrixXf GRM = 1/ m * genos.transpose() * genos;
    Eigen::VectorXf diag = GRM.diagonal();
    float mean_diag = diag.mean();
    cerr << "Rescaling by mean of GRM diagonal: " << mean_diag << endl;
    GRM = 1/mean_diag * GRM;
    return GRM;
}

void write_pheno_GRM(
    const Eigen::MatrixXf &GRM,
    const vector<int> &pheno,
    const char *grm_out,
    const char *pheno_out
) {
    cerr << "Writing phenotypes to " << pheno_out << endl;
    std::ofstream pheno_out_stream;
    pheno_out_stream.open(pheno_out); assert(pheno_out_stream);
    for(int i = 0; i < pheno.size(); i++) {
        pheno_out_stream << i+1 << "\t" << i+1 << "\t" << pheno[i] << endl;
    }
    pheno_out_stream.close();

    string grm_id_out = string(grm_out) + ".grm.id";
    cerr << "Writing GRM ID file to " << grm_id_out << endl;
    std::ofstream grm_id_out_stream;
    grm_id_out_stream.open(grm_id_out); assert(grm_id_out_stream);
    for (int i = 0; i < pheno.size(); i++) {
        grm_id_out_stream << i+1 << "\t" << i+1 << endl;
    }
    grm_id_out_stream.close();

    // Write GRM binary files
    string grm_bin_out = string(grm_out) + ".grm.bin";
    cerr << "Writing GRM binary files to " << grm_bin_out << endl;

    std::ofstream grm_bin_stream(grm_bin_out, std::ios::binary);
    assert(grm_bin_stream);
    for (int i = 0; i < GRM.rows(); i++) {
        for (int j = 0; j <= i; j++) {
            float grm_value = GRM(i, j);
            grm_bin_stream.write(reinterpret_cast<const char*>(&grm_value), sizeof(float));
        }
    }
    grm_bin_stream.close();
}

Eigen::MatrixXf read_difference_genotypes(
    htsFile *fin, 
    bcf_hdr_t *hdr, 
    vector<queue<CNV>> &mCALookupArray, 
    const vector<bool> &removeArray, 
    const set<string> &discardSet,
    vector<int> &pheno
) {
    bcf1_t *rec = bcf_init();
    int *GT = NULL, GT_buf_size = 0;
    int *AS = NULL, AS_buf_size = 0;
    int chr;
    int progress = 0;
    time_t timestamp;
    char *rid_to_chr = mapContigToChrom(hdr);
    int n = bcf_hdr_nsamples(hdr);
    vector<int> mCA_phase(n, bcf_int32_missing);
    vector<bool> keepArray(n, false);
    vector<vector<float>> difference_genos;
    int n_missing = 0;
    int m = 0;

    while (bcf_read(fin, hdr, rec) == 0) {
        progress += 1;
        chr = rid_to_chr[rec->rid];
        if (chr == 0) continue;
        assert(rec->n_allele == 2);

        int nGT = bcf_get_genotypes(hdr, rec, &GT, &GT_buf_size);
        if (nGT != 2 * n) {
            // cerr << "WARNING: missing genotypes; skipping chr" << chr << ":" << rec->pos + 1 << endl;
            continue;
        }
        assert(bcf_unpack(rec, BCF_UN_STR) == 0);
        char var_name_buffer[256];
        if (chr != 23 ) {
            snprintf(var_name_buffer, sizeof(var_name_buffer), "%d:%ld:%s:%s", chr, rec->pos + 1, rec->d.allele[0], rec->d.allele[1]);
        } else {
            snprintf(var_name_buffer, sizeof(var_name_buffer), "X:%ld:%s:%s", rec->pos + 1, rec->d.allele[0], rec->d.allele[1]);
        }
        string var_name(var_name_buffer);
        if (discardSet.count(var_name) == 1) continue;

        int nAS = bcf_get_format_int32(hdr, rec, "AS", &AS, &AS_buf_size);

        if (AS_buf_size == 0) {
            AS = (int *)malloc(n * sizeof(int));
            AS_buf_size = n;
        }
        if (nAS != n) {
            nAS = n;
            for (int i = 0; i < n; i++) {
                if (mCA_phase[i] == bcf_int32_missing || !variantInCNV(chr, rec->pos + 1, mCALookupArray[i])) 
                    AS[i] = 0;
                else 
                    AS[i] = ALLELIC_SHIFT_IS_1(mCA_phase[i], GT[2 * i]) ? 1 : -1;
            }
        } else {
            for (int i = 0; i < n; i++) {
                if (AS[i] != -1 && AS[i] != 1) 
                    AS[i] = variantInCNV(chr, rec->pos + 1, mCALookupArray[i]) 
                            ? (ALLELIC_SHIFT_IS_1(mCA_phase[i], GT[2 * i]) ? 1 : -1) 
                            : bcf_int32_missing;
                else 
                    mCA_phase[i] = MCA_PHASE_IS_1(AS[i], GT[2 * i]) ? 1 : 0;
            }
        }

        float missing = 0, AN = 0, AC = 0;
        for (int i = 0; i < n; i++) {
            if (removeArray[i]) continue;
            if (bcf_gt_is_missing(GT[2 * i]) || bcf_gt_is_missing(GT[2 * i + 1])) {
                missing += 2;
            } else {
                AC += bcf_gt_allele(GT[2 * i]) + bcf_gt_allele(GT[2 * i + 1]);
                AN += 2;
            }
        }
        if(AN == 0 || AC == 0 || AC == AN) continue;
        if ((missing / (missing + AN)) > 0.05) {
            // cerr << "WARNING: missing " << missing/(missing+AN) << " ; skipping chr" << chr << ":" << rec->pos + 1 << endl;
            n_missing++;
            continue;
        }
        float AF = AC / AN;

        vector<float> variant_difference_genos(n, 0);
        bool keep_variant = false;
        for (int i = 0; i < n; i++) {
            int offset = i % 2;
            if (removeArray[i]) continue;
            if (AS[i] != 1 && AS[i] != -1) continue;
            if (!variantInCNV(chr, rec->pos + 1, mCALookupArray[i])) continue;
            if (bcf_gt_is_missing(GT[2 * i]) || bcf_gt_is_missing(GT[2 * i + 1])) continue;
            assert(AS[i] == 1 || AS[i] == -1);
            float diff  = MCA_PHASE_IS_1(AS[i], GT[2 * i + offset]) 
                         ? bcf_gt_allele(GT[2 * i + 1]) - bcf_gt_allele(GT[2 * i]) 
                         : bcf_gt_allele(GT[2 * i]) - bcf_gt_allele(GT[2 * i + 1]);
            assert(diff == 1 || diff == 0 || diff == -1);
            assert(AF > 0 && AF < 1);
            diff /= sqrt(2 * AF * (1 - AF));
            variant_difference_genos[i] = diff;
            pheno[i] = offset;
            keepArray[i] = true;
            keep_variant = true;
        }
        if(keep_variant) {
            difference_genos.push_back(variant_difference_genos);
            m += 1;
            if (m % 1000 == 0) {
                time(&timestamp);
                cerr << strtok(ctime(&timestamp), "\n") << ": chr" << chr << ":" << rec->pos+1 ;
                cerr << " (" << m << " variants kept; ";
                cerr <<  progress << " variants skipped)" << endl;
            }
        }

    }
    cerr << "Skipped " << n_missing << " variants due to missingness" << endl;
    assert(difference_genos.size() > 0);
    int num_ind = 0;
    int skipped_ind = 0;
    vector<int> ind_map(n, -1);
    for(int i = 0; i < n; i++) {
        if(keepArray[i]) num_ind++;
        else {
            skipped_ind++;
            continue;
        }
        ind_map[i] = i - skipped_ind;
        pheno[ind_map[i]] = pheno[i];
    }
    assert(num_ind > 0);
    pheno.resize(num_ind);

    // Convert difference_genos to an Eigen matrix
    int num_var = difference_genos.size();
    cerr << "Read " << num_var << " variants for " << num_ind << " individuals" << endl;
    Eigen::MatrixXf geno_matrix(num_var, num_ind);
    for (int i = 0; i < num_var; i++) {
        for (int j = 0; j < n; j++) {
            if (!keepArray[j]) continue;
            assert(ind_map[j] != -1);
            assert(ind_map[j] <= j);
            assert(ind_map[j] < num_ind);
            geno_matrix(i, ind_map[j]) = difference_genos[i][j];
        }
    }

    free(AS);
    free(rid_to_chr);
    bcf_destroy(rec);

    return geno_matrix;
}

int write_CNLOH_GRM_main(std::vector<std::string> &argv) {
    const char *calls_fin = NULL;
    const char *remove_fin = NULL;
    const char *grm_out = NULL;
    const char *pheno_out = NULL;
    const char *discard_fin = NULL;
    string target_chrom_arm;
    string event_type = "CN-LOH";

    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "produce help message")
        ("calls", po::value<std::string>(), "file with the mCA calls")
        ("grm-out", po::value<std::string>(), "file to write the GRM to")
        ("pheno-out", po::value<std::string>(), "file to write the phenotype to")
        ("chrom-arm", po::value<std::string>(), "chromosome arm to analyze")
        ("discard", po::value<std::string>(), "file with pruned variants to discard from the analysis")
        ("cn", po::value<std::string>(), "event type to test")
        ("remove", po::value<std::string>(), "file with individuals to remove from the analysis");

    po::variables_map vm;
    try {
        po::store(po::command_line_parser(argv).options(desc).run(), vm);
        po::notify(vm);
    } catch(exception& e) {
        cerr << "error: " << e.what() << endl;
        return 1;
    } catch(...) {
        cerr << "Exception of unknown type!" << endl;
        return 1;
    }

    if (vm.count("help")) {
        cerr << desc << endl;
        return 0;
    }
    if (vm.count("calls")) {
        cerr << "Reading calls from " << vm["calls"].as<std::string>() << endl;
        calls_fin = vm["calls"].as<std::string>().c_str();
    } else {
        cerr << "Error: need to supply a file with calls" << endl;
        cerr << desc << endl;
        return 1;
    }	
    if (vm.count("chrom-arm")) {
        target_chrom_arm = vm["chrom-arm"].as<std::string>();
    } 
    else {
        cerr << "Error: need to supply a chromosome arm" << endl;
        cerr << desc << endl;
        return 1;
    }
    if (vm.count("grm-out")) {
        grm_out = vm["grm-out"].as<std::string>().c_str();
    } else {
        grm_out = "GRM.grm";
    }	
    if (vm.count("pheno-out")) {
        pheno_out = vm["pheno-out"].as<std::string>().c_str();
    } else {
        pheno_out = "pheno.txt";
    }	
    if(vm.count("discard")) {
        cerr << "Discarding variants in " << vm["discard"].as<std::string>() << endl;
        discard_fin = vm["discard"].as<std::string>().c_str();
    }
    if (vm.count("remove")) {
        cerr << "Removing individuals in " << vm["remove"].as<std::string>() << endl;
        remove_fin = vm["remove"].as<std::string>().c_str();
    }
    if (vm.count("cn")) {
        event_type = vm["cn"].as<std::string>();
    }
    cerr << "Testing " << event_type << endl;
    if(event_type != "CN-LOH" && event_type != "GAIN" && event_type != "LOSS") {
        cerr << "Error: must specify copy-number as either CN-LOH, GAIN, or LOSS" << endl;
        return 1;
    }

    htsFile *fin = hts_open("-", "r");
    assert(fin != NULL);
    bcf_hdr_t *hdr = bcf_hdr_read(fin);
    assert(hdr != NULL);
    int n = bcf_hdr_nsamples(hdr);

    vector<bool> removeArray(n, false);
    readRemoveArray(remove_fin, hdr, removeArray);

    set<string> discardSet;
    readDiscardSet(discard_fin, discardSet);

    vector<queue<CNV>> mCALookupArray(n);
    int num_ind_with_calls = read_mCAFile(calls_fin, event_type, target_chrom_arm, removeArray, hdr, mCALookupArray);
    if(num_ind_with_calls == 0) {
        cerr << "No individuals with " << event_type << " calls on " << target_chrom_arm << endl;
        return 0;
    }

    vector<int> mCA_phase(n, bcf_int32_missing);

    vector<int> pheno(n, -1);
    cerr << "Reading genotypes..." << endl;
    Eigen::MatrixXf difference_genos = read_difference_genotypes(fin, hdr, mCALookupArray, removeArray, discardSet, pheno);
    bcf_hdr_destroy(hdr);
    hts_close(fin);

    cerr << "Computing GRM..." << endl;
    Eigen::MatrixXf GRM = compute_GRM(difference_genos);
    write_pheno_GRM(GRM, pheno, grm_out, pheno_out);

    return 0;
}

int main(int argc, char **argv) {
    std::vector<std::string> args;
    for (int i = 1; i < argc; ++i) args.emplace_back(argv[i]);
    return write_CNLOH_GRM_main(args);
}

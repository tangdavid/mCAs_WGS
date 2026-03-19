#include <htslib/vcf.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <boost/program_options.hpp>
#include <boost/math/distributions/binomial.hpp>
#include "prepare_bcf.h"
#include <htslib/kfunc.h>
#include <ctime>
#include <set>
#include "gwas_tools.h"

using namespace std;
namespace po=boost::program_options;


// mCA_phase[i] == 1 if AS[i] == -1 and GT[2*i] == 0 or AS[i] == 1 and GT[2*i] == 1; else 0
// AS[i] == -1 means that mCA_phase == 0 if GT[2*i] == 1 and mCA_phase == 1 if GT[2*i] == 0; else 1
#define MCA_PHASE_IS_1(as,first_allele) ((as == -1 && bcf_gt_allele(first_allele) == 0) || (as == 1 && bcf_gt_allele(first_allele) == 1))
#define ALLELIC_SHIFT_IS_1(phase,first_allele) ((phase == 0 && bcf_gt_allele(first_allele) == 0) || (phase == 1 && bcf_gt_allele(first_allele) == 1))

int compute_cis_fisher_p_main(std::vector<std::string> &argv) {

    const char *calls_fin = NULL;
    const char *remove_fin = NULL;
    string event_type = "CN-LOH";

    po::variables_map vm;
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "produce help message")
        ("calls", po::value<std::string>(), "file with the mCA calls\n")
        ("remove", po::value<std::string>(), "file with individuals to remove from the analysis\n")
        ("cn", po::value<std::string>(), "event type to test\n");

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
    if (vm.count("calls")) {
        cerr << "Reading calls from " << vm["calls"].as<std::string>() << endl;
        calls_fin = vm["calls"].as<std::string>().c_str();
    } else {
        cerr << "Error: need to supply a file with calls" << endl;
        return 1;
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
    htsFile *fin = hts_open("-", "r"); assert(fin != NULL);
    bcf_hdr_t *hdr = bcf_hdr_read(fin); assert(hdr != NULL);
    int n = bcf_hdr_nsamples(hdr);

    char *rid_to_chr = mapContigToChrom(hdr);
    bcf1_t *rec = bcf_init();

    vector<bool> removeArray(n, false);
    int n_removed = readRemoveArray(remove_fin, hdr, removeArray);

    vector<queue<CNV>> mCALookupArray(n);
    string target_chrom_arm = "";
    read_mCAFile(calls_fin, event_type, target_chrom_arm, removeArray, hdr, mCALookupArray);

    int *GT = NULL, GT_buf_size=0;
    int *AS = NULL, AS_buf_size=0;
    float *PP = NULL; int PP_buf_size=0;
    int chr;
    int progress = 0;
    time_t timestamp;
    vector<int> mCA_phase(n, bcf_int32_missing);

    cout << "chr\tpos\tvariant\tref\talt\tcarrier_mCA\tnoncarrier_mCA\tcarrier_control\tnoncarrier_control\tfisher_p\tcarrier_overrep\tcarrier_underrep\tbad_phase\tbinom_p\twarnings" << endl;

    while (bcf_read(fin, hdr, rec) == 0) {
        chr = rid_to_chr[rec->rid];
        if(chr == 0) continue;
        assert(rec->n_allele == 2);
        int nGT = bcf_get_genotypes(hdr, rec, &GT, &GT_buf_size);
        if(nGT != 2*n) {
            cerr << "WARNING: missing genotypes; skipping chr" << chr << ":" << rec->pos+1 << endl;
            continue;
        }
        assert(bcf_unpack(rec, BCF_UN_STR) == 0);

        time(&timestamp);
        if (progress % 1000 == 0) cerr << strtok(ctime(&timestamp), "\n") << ": chr" << chr << ":" << rec->pos+1 << " (" << progress << " variants processed)" << endl;
        progress += 1;

        int nAS = bcf_get_format_int32(hdr, rec, "AS", &AS, &AS_buf_size);

        if (AS_buf_size == 0) {
            AS = (int*)malloc(n*sizeof(int));
            AS_buf_size = n;
        }   
        if (nAS != n) {
            // set AS field based on mCA_phase and GT
            nAS = n;
            for (int i = 0; i < n; i++) {
                if (mCA_phase[i] == bcf_int32_missing || !variantInCNV(chr, rec->pos+1, mCALookupArray[i])) AS[i] = 0;
                else AS[i] = ALLELIC_SHIFT_IS_1(mCA_phase[i], GT[2*i]) ? 1 : -1;
            }
            
        }  else {
            // update the mCA_phase field based on AS and GT
            for (int i = 0; i < n; i++) {
                if(AS[i] != -1 && AS[i] != 1) AS[i] = variantInCNV(chr, rec->pos+1, mCALookupArray[i]) ? (ALLELIC_SHIFT_IS_1(mCA_phase[i], GT[2*i]) ? 1 : -1) : bcf_int32_missing;
                else mCA_phase[i] = MCA_PHASE_IS_1(AS[i], GT[2*i]) ? 1 : 0;
            }
        }
        
        int nPP = bcf_get_format_float(hdr, rec, "PP", &PP, &PP_buf_size);       
        //if(nPP != n) continue;

        int c_mCA = 0, nc_mCA = 0, c_control = 0, nc_control = 0;
        int bad_phase = 0, overrepresented = 0, underrepresented = 0;
        
        for (int i=0; i < n; i++) {
            if(removeArray[i]) continue;
            int is_het = (bcf_gt_allele(GT[2*i]) + bcf_gt_allele(GT[2*i+1]) == 1) ? 1: 0;
            if(variantInCNV(chr, rec->pos+1, mCALookupArray[i])) {
                nc_mCA += 1 - is_het;
                c_mCA += is_het;
                if (is_het != 1) continue;
                if (nPP == n && PP[i] < 0.8) bad_phase += 1;
                else {
                    // AS < 0 means allelic shift in direction of the alt variant
                    overrepresented += (AS[i] < 0) ? 1 : 0;
                    underrepresented += (AS[i] > 0) ? 1 : 0;
                }
            } else {
                nc_control += 1 - is_het;
                c_control += is_het;
            }
        }
        //assert(bad_phase + overrepresented + underrepresented == c_mCA);
        string warning = ".";
        if(bad_phase + overrepresented + underrepresented != c_mCA) 
            warning="missing_skew";
        assert(c_mCA + nc_mCA + c_control + nc_control == n - n_removed);

        double left_p, right_p, two_p;
        kt_fisher_exact(c_mCA, nc_mCA, c_control, nc_control, &left_p, &right_p, &two_p);

        boost::math::binomial_distribution<> binom(overrepresented+underrepresented,0.5);
        double binom_p = boost::math::cdf(binom, min(underrepresented, overrepresented))*2;
        binom_p = min(1.0, binom_p);
        
        cout << "chr" << chr << "\t" << rec->pos+1 << "\t" << rec->d.id << "\t" << rec->d.allele[0] << "\t" << rec->d.allele[1] << "\t"
        << c_mCA << "\t" << nc_mCA << "\t" << c_control << "\t" << nc_control << "\t" << right_p << "\t"
        << overrepresented << "\t" << underrepresented << "\t" << bad_phase << "\t" << binom_p << "\t" << warning << endl;
    }
    free(GT); free(AS); free(PP);
    free(rid_to_chr);

    bcf_destroy(rec);
    bcf_hdr_destroy(hdr);
    hts_close(fin);

    return 0;
}

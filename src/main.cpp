#include <iostream>
#include <cstring>
#include <vector>

using namespace std;

int prepare_bcf_main(std::vector<std::string> &argv);
int make_calls_main(std::vector<std::string> &argv);
int compute_cis_fisher_p_main(std::vector<std::string> &argv);
int compute_polygenic_drive_main(std::vector<std::string> &argv);

void printModes(){
    cerr << "Usage:" << endl;
    cerr << "  mCAs_WGS [mode] [options]" << endl;
    cerr << "  eg: mCAs_WGS make-calls --help" << endl;
    cerr << "Available modes:" << endl;
    cerr << "  prepare-bcf\t\t\tTools to merge in ADs, mask CNVs, perform genotype QC, etc." << endl;
    cerr << "  make-calls\t\t\tCall mCAs based on BAF" << endl;
    cerr << "  compute-cis-fisher\t\tPerform GWAS for cis allelic shift" << endl;
    cerr << "  compute-polygenic-drive\tCompute differential PRS within CN-LOH boundaries" << endl;
}

int main(int argc, char *argv[]) {
	std::vector < std::string > args;
    if (argc < 2){
        printModes();
        exit(0);
    }
	for (int a = 2 ; a < argc ; a ++) args.push_back(std::string(argv[a]));

    if (strcmp(argv[1], "prepare-bcf") == 0) prepare_bcf_main(args);
    else if (strcmp(argv[1], "make-calls") == 0) make_calls_main(args);
    else if (strcmp(argv[1], "compute-cis-fisher") == 0) compute_cis_fisher_p_main(args);
    else if (strcmp(argv[1], "compute-polygenic-drive") == 0) compute_polygenic_drive_main(args);
    else if (strcmp(argv[1], "--help") == 0) {
        printModes();
        exit(0);
    } else {
        printModes();
        cerr << "Error: Unrecognized mode " << argv[1] << endl;
        exit(1);
    }
}

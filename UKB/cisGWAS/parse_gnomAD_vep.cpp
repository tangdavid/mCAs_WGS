#include <iostream>
#include <cstdio>
#include <cstring>
#include <cstdlib>

using namespace std;

char vep[1000000];

int main(void) {
    string chrom, pos, id, ref, alt, qual, filter; int AN, AN_non_ukb, AN_non_ukb_nfe, AC, AC_non_ukb, AC_non_ukb_nfe;
    while (cin >> chrom >> pos >> id >> ref >> alt >> qual >> filter
         >> AN >> AN_non_ukb >> AN_non_ukb_nfe >> AC >> AC_non_ukb >> AC_non_ukb_nfe) {
        scanf("%s", vep);
        if (filter != "PASS") continue;
        int i = 0;
        bool good = true;
        bool is_mane_select = false;
        bool is_LoF = false;
        string Consequence, SYMBOL, IMPACT, HGVSp/*, Gene, Feature, LoF, LoF_filter, LoF_flags*/;
        char *token, *str, *tofree;
        tofree = str = strdup(vep);
        while ((token = strsep(&str, "|")) != NULL) {
            i++;
            switch (i) {
                case 2: Consequence = token; if (strlen(token)==0) good = false; break;
                case 3: if (strcmp(token, "HIGH")!=0 && strcmp(token, "MODERATE")!=0) good = false; break; // IMPACT
                case 4: SYMBOL = token; if (strlen(token)==0) good = false; break;
                case 6: if (strcmp(token, "Transcript")!=0) good = false; break; // Feature_type
                case 8: if (strcmp(token, "protein_coding")!=0) good = false; break; // BIOTYPE
                case 12: HGVSp=token; break;
                case 26: if (strlen(token)!=0) is_mane_select = true; break; // MANE_SELECT
                case 33: if (strcmp(token, "Ensembl")!=0) good = false; break; // SOURCE
                case 43: if (strcmp(token, "HC")==0) is_LoF = true; break; // LoF
                case 44: if (strlen(token)!=0) is_LoF = false; break; // LoF_filter
                case 45: if (strlen(token)!=0) is_LoF = false; break; // LoF_flags
            }
            if (i == 45) { // reached LoF_flags, which merges with the next vep record
                if (good) {
                    HGVSp = HGVSp.size() == 0 ? "." : HGVSp;
        
                    cout << chrom << "\t" << pos << "\t" << id << "\t" << ref << "\t" << alt << "\t" 
                      << AN-AN_non_ukb << "\t" << AC-AC_non_ukb << "\t" << AN_non_ukb_nfe << "\t" << AC_non_ukb_nfe << "\t" 
                      << SYMBOL << "\t" << Consequence << "\t" << HGVSp << "\t" 
                      << (is_mane_select ? "MANE_SELECT" : ".") << "\t" << (is_LoF ? "LoF" : ".") << endl;
                }
                i = 0; // reset counter
                good = true;
                is_mane_select = false;
                is_LoF = false;
            }
        }
        free(tofree);
    }
    return 0;
}

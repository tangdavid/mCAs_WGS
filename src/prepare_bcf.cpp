#include "prepare_bcf.h"
#include <boost/program_options.hpp>

using namespace std;
namespace po=boost::program_options;

int chr2int(string chr) {
    if (chr == "chrX") {
        return 23;
    }
    chr.erase(0, 3);
    int res;
    try {
        res = stoi(chr);
        return res;
    }
    catch (exception &err) {
        return -1;
    }
}

bool variantInCNV(
    int chrNum, 
    uint32_t pos, 
    queue<CNV> &q
) {
    bool res = false;
    CNV currentCNV;

    while(!q.empty()) {
        currentCNV = q.front();
        if (currentCNV.chrNum > chrNum) {
            break;
        }
        if (currentCNV.chrNum < chrNum) {
            q.pop();
            continue;
        }
        if (pos > currentCNV.end) {
            q.pop();
            continue;
        }
        res = pos >= currentCNV.start;
        break;
    }

    return res;
}

void maskCNV(
    int chrNum, 
    uint32_t pos, 
    const bcf_hdr_t *hdr, 
    vector<queue<CNV>> &table, 
    int *fmt_field
) {
    for (int i = 0; i < bcf_hdr_nsamples(hdr); i++) {
        if (variantInCNV(chrNum, pos, table[i])) {
            fmt_field[2*i] = fmt_field[2*i+1] = bcf_int32_missing;
        }
    }
}

void adjustRefBias(
    int n, 
    double bias, 
    int *fmt_AD
) {
    for (int i = 0; i < n; i++) {
        if (
            abs(bias-0.5) >= 0.05 || 
            fmt_AD[2*i] == bcf_int32_missing || 
            fmt_AD[2*i+1] == bcf_int32_missing
        ) {
            fmt_AD[2*i] = fmt_AD[2*i+1] = bcf_int32_missing;
        } 
        else {
            if (bias > 0.5) {
                fmt_AD[2*i] = round(fmt_AD[2*i] * (1-bias)/bias);
            }
            else {
                fmt_AD[2*i+1] = round(fmt_AD[2*i+1] * bias/(1-bias));
            }
        }
    }
}

void computeRefBias(const int *ADs, int n, const bcf1_t *rec, ofstream *file) {
    int numSample = 0;
    int refCounts = 0;
    int altCounts = 0;
    for (int i = 0; i < n; i++) {
        if (ADs[2*i] == bcf_int32_missing || ADs[2*i+1] == bcf_int32_missing) {
            continue;
        }
        refCounts += ADs[2*i];
        altCounts += ADs[2*i+1];
        numSample += 1;
    }
    if (numSample == 0) {
        return;
    }
    *file << "chr" << rec->rid+1 << "\t"; 
    *file << rec->pos+1 << "\t"; 
    *file << rec->d.allele[0] << "\t"; 
    *file << rec->d.allele[1] << "\t"; 
    *file << refCounts << "\t"; 
    *file << altCounts << "\t"; 
    *file << numSample; 
    *file << endl; 
}

map<string, double> parseRefBias(ifstream *file) {
    string chr_str;
    uint32_t pos;
    string ref;
    string alt;
    double ref_count;
    double alt_count;
    char variant[1000];
    int max_len =sizeof(variant);

    map<string, double> table; 

    while(*file >> chr_str >> pos >> ref >> alt >> ref_count >> alt_count) {
        int chrNum = chr2int(chr_str);
        int b = snprintf(variant, max_len, "%d_%d_%s_%s", chrNum, pos, ref.c_str(), alt.c_str());
        if (b >= max_len) {
            cerr << "WARNING: variant name buffer length exceeded when parsing ref bias" << endl;
        }
        table[variant] = (ref_count+alt_count == 0) ? -1 : ref_count/(ref_count + alt_count);
    }
    return table;
}

map<string, queue<CNV>> parseCNVs(ifstream *file, bool mask_deletions) {
    string sample;
    string chr_str;
    uint32_t start;
    uint32_t end;
    string type;

    map<string, queue<CNV>> table;

    while(*file >> sample >> chr_str >> start >> end >> type) {
        if (type == "DEL" && !mask_deletions) {
            continue;
        }   
        int chrNum = chr2int(chr_str);
        if (chrNum < 0) {
            continue;
        }
        CNV record = {chrNum, start, end};
        if (!table[sample].empty()) {
            const CNV &last_record = table[sample].back();
            assert(last_record.chrNum <= chrNum && "CNVs are not sorted by chromosome");
            if(last_record.chrNum == chrNum) {
                assert(last_record.start <= record.start && "CNVs are not sorted by start position");
            }
        }
        table[sample].push(record);
    }
    return table;
}

variantAD_t consumeVariants(ifstream *f, const bcf1_t *rec) {
    variantAD_t var;
    var.pos = 1e9;
    var.rid = 99;
    string chr_str;
    
    // discard variants that are not in the phased file 
    while(*f >> chr_str >> var.pos >> var.ref >> var.alt >> var.ad1 >> var.ad2) {
        var.rid = chr2int(chr_str)-1;
        var.pos -= 1;
        if (var.rid > rec->rid || (var.rid == rec->rid && var.pos >= rec->pos)) break;
    }
    return var;
}

int checkVariantEquality(const variantAD_t &var_AD, const bcf1_t* var_BCF) {
    return (
        (var_AD.rid == var_BCF->rid) &&
        (var_AD.pos == var_BCF->pos) && 
        (var_AD.ref == var_BCF->d.allele[0]) && 
        (var_AD.alt == var_BCF->d.allele[1])
    );
}

int writeADs(int i, int *AD, const variantAD_t &var_AD, const int *GT) {
    int gt1 = bcf_gt_allele(GT[2*i]);
    int gt2 = bcf_gt_allele(GT[2*i+1]);
    
    if ((gt1 + gt2) != 1) {
        return 1; // keep track of how many variants are discarded
    }
    
    AD[2*i] = var_AD.ad1;
    AD[2*i + 1] = var_AD.ad2;
    return 0;
}

int parseChrStrToChr0(const char *chrStr) {
    int chr0 = 0;
    if (strcmp(chrStr, "chrX")==0)
        chr0 = 23-1;
    else if (strcmp(chrStr, "chrY")==0)
        chr0 = 24-1;
    else {
        assert(strncmp(chrStr, "chr", 3)==0);
        sscanf(chrStr, "chr%d", &chr0); chr0--;
        assert(chr0 >= 0 && chr0 < 22);
    }
    return chr0;
}


char *mapContigToChrom(const bcf_hdr_t *hdr) {
    char *rid_to_chr = (char *) calloc(hdr->n[BCF_DT_CTG], 1);
    for (int rid = 0; rid < hdr->n[BCF_DT_CTG]; rid++) {
        const char *contig_name = hdr->id[BCF_DT_CTG][rid].key;
        if(strlen(contig_name) > 3 && contig_name[0]=='c' && contig_name[1]=='h' && contig_name[2]=='r') {
            if(strlen(contig_name)==4) {
                if (isdigit(contig_name[3])) 
                    rid_to_chr[rid] = contig_name[3]-'0';
                else if (contig_name[3]=='X') 
                    rid_to_chr[rid] = 23;
                else if (contig_name[3]=='Y')
                    rid_to_chr[rid] = 24;
            }
            else if (strlen(contig_name) == 5 && isdigit(contig_name[3]) && isdigit(contig_name[4])) {
                rid_to_chr[rid] = 10*(contig_name[3]-'0') + contig_name[4]-'0';
            }
        }
    }
    return rid_to_chr;
}

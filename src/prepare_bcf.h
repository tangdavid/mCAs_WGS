#ifndef __PREPARE_BCF_H__
#define __PREPARE_BCF_H__

#include <htslib/vcf.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <queue>
#include <map>
#include <getopt.h>
#include <math.h>
#include <string>

using namespace std;

struct CNV {
    int chrNum;
    uint32_t start;
    uint32_t end;
};

struct variantAD_t {
    int rid;
    int pos;
    string ref;
    string alt;
    int ad1;
    int ad2;
} ;

variantAD_t consumeVariants(ifstream *f, const bcf1_t* rec);

int checkVariantEquality(const variantAD_t &var_AD, const bcf1_t* var_BCF);

int writeADs(int i, int *AD, const variantAD_t &var_AD, const int *GT);

/**
 * @brief checks whether or not a variant is in a CNV and updates
 * expects CNVs to be stored in order in a queue 
 * 
 * @param rid is chromosome number
 * @param pos is variant position
 * @param q is queue of all CNVs called in the sample
 * @return true if variant is in a CNV
 * @return false otherwise
 */
bool variantInCNV(int rid, uint32_t pos, queue<CNV> &q);

/**
 * @brief for a given variant, masks fmt_fields of samples with a CNV 
 * overlapping the variant position
 * 
 * @param rid is chromosome number
 * @param pos is variant position
 * @param hdr is the bcf header
 * @param table is list {queue of all CNVs called in the sample}
 * @param fmt_field is a list of ADs at the given position
 */
void maskCNV(int rid,uint32_t pos, const bcf_hdr_t *hdr, vector<queue<CNV>> &table, int *fmt_field);

/**
 * @brief convert chromosome name to number; chrX -> is considered 23
 * 
 * @param chr is chromosome name
 * @return int chromosome number
 */
int chr2int(string chr);

void adjustRefBias(int n, double bias, int *fmt_AD);

void computeRefBias(const int *ADs, int n, const bcf1_t *rec, ofstream *file);

map<string, double> parseRefBias(ifstream *file);

map<string, queue<CNV>> parseCNVs(ifstream *file, bool mask_deletions=false);

int parseChrStrToChr0(const char *chrStr);

char *mapContigToChrom(const bcf_hdr_t *hdr);

#endif

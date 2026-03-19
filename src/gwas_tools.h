#ifndef GWAS_TOOLS_H
#define GWAS_TOOLS_H

#include <string>
#include <vector>
#include <queue>
#include <map>
#include <set>
#include <cassert>
#include <fstream>
#include <sstream>
#include <iostream>
#include "prepare_bcf.h"

int read_mCAFile(
    const char *calls_fin,
    const std::string &event_type,
    const std::string &target_chrom_arm, // Optional: empty string means no filtering by chromosome arm
    const std::vector<bool> &removeArray,
    const bcf_hdr_t *hdr,
    std::vector<std::queue<CNV>> &mCALookupArray
);

int readRemoveArray(
    const char *remove_fin,
    const bcf_hdr_t *hdr,
    std::vector<bool> &removeArray
);

int readDiscardSet(
    const char *discard_fin,
    std::set<std::string> &discardSet
);

#endif // GWAS_TOOLS_H

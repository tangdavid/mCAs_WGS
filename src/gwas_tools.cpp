#include "gwas_tools.h"
#include "genome_rules.h"

int read_mCAFile(
    const char *calls_fin,
    const std::string &event_type,
    const std::string &target_chrom_arm,
    const std::vector<bool> &removeArray,
    const bcf_hdr_t *hdr,
    std::vector<std::queue<CNV>> &mCALookupArray
) {
    std::ifstream mCAFile(calls_fin);
    assert(mCAFile);
    std::map<std::string, std::queue<CNV>> mCALookup;
    std::string header_line;
    getline(mCAFile, header_line); // Read the header line
    std::vector<std::string> headers;
    std::stringstream header_stream(header_line);
    std::string column;
    while (getline(header_stream, column, '\t')) {
        headers.push_back(column); // Store column names
    }

    genome_rules_t *genome_rules;
    char build[] = "GRCh38";
    genome_rules = genome_init_alias(stderr, build, hdr);
    readlist_short_arms(genome_rules, "13,14,15,21,22,chr13,chr14,chr15,chr21,chr22", hdr);

    std::map<std::string, int> chrToContig;
    for (int rid = 0; rid < hdr->n[BCF_DT_CTG]; rid++) {
        string chr_str(hdr->id[BCF_DT_CTG][rid].key);
        chrToContig[chr_str]=rid;
    }

    // Identify relevant column indices
    int sample_col = -1, chr_col = -1, start_col = -1, end_col = -1, type_col = -1, p_col = -1, q_col = -1;
    for (size_t i = 0; i < headers.size(); i++) {
        if (headers[i] == "ID") sample_col = i;
        else if (headers[i] == "chr") chr_col = i;
        else if (headers[i] == "bpStart") start_col = i;
        else if (headers[i] == "bpEnd") end_col = i;
        else if (headers[i] == "type") type_col = i;
        else if (headers[i] == "p") p_col = i;
        else if (headers[i] == "q") q_col = i;
    }
    assert(sample_col != -1 && chr_col != -1 && start_col != -1 && end_col != -1 && type_col != -1);

    std::string line;
    while (getline(mCAFile, line)) {
        std::stringstream line_stream(line);
        std::vector<std::string> fields;
        while (getline(line_stream, column, '\t')) {
            fields.push_back(column); // Split line into fields
        }
        std::string sample = fields[sample_col];
        std::string chr_str = fields[chr_col];
        if(!chrToContig.count(chr_str)) continue;

        int start = (fields[p_col] == "T") ? 1 : std::stoi(fields[start_col]);
        int end = (fields[q_col] == "T") ? 1000000000 : std::stoi(fields[end_col]);
        start = (abs(start - genome_rules->cen_end[chrToContig[chr_str]])<500000) ? genome_rules->cen_end[chrToContig[chr_str]] : start;
        end = (abs(end - genome_rules->cen_beg[chrToContig[chr_str]])<500000) ? genome_rules->cen_beg[chrToContig[chr_str]] : end;
        assert(start > 0 && end > 0);

        std::string type = fields[type_col];
        std::string chr_arm = (fields[p_col] != "N") ? ((fields[q_col] != "N") ? "pq" : "p") : "q";

        if (type != event_type) continue;
        if (!target_chrom_arm.empty() && chr_str + chr_arm != target_chrom_arm) continue;
        int chrNum = chr2int(chr_str);
        CNV record = {chrNum, static_cast<uint32_t>(start), static_cast<uint32_t>(end)};
        if (!mCALookup[sample].empty()) {
            const CNV &last_record = mCALookup[sample].back();
            assert(last_record.chrNum <= chrNum && "calls are not sorted by chromosome");
            if(last_record.chrNum == chrNum) {
                assert(last_record.start <= record.start && "calls are not sorted by start position");
            }
        }
        mCALookup[sample].push(record);
    }

    int n = bcf_hdr_nsamples(hdr);
    int num_ind = 0;
    for (int i = 0; i < n; i++) {
        if (removeArray[i]) continue;
        mCALookupArray[i] = mCALookup[hdr->samples[i]];
        if (mCALookupArray[i].size() > 0) num_ind++;
    }
    std::cerr << "Using " << num_ind << " individuals with " << event_type << " calls";
    if (!target_chrom_arm.empty()) {
        std::cerr << " on " << target_chrom_arm;
    }
    std::cerr << std::endl;
    mCAFile.close();
    return num_ind; 
}

int readRemoveArray(
    const char *remove_fin,
    const bcf_hdr_t *hdr,
    std::vector<bool> &removeArray
) {
    if (!remove_fin) return 0;

    std::ifstream removeFile(remove_fin);
    assert(removeFile);
    std::set<std::string> removeSet;
    std::string sample;

    while (removeFile >> sample) {
        removeSet.insert(sample);
    }

    int n_removed = 0;
    for (int i = 0; i < removeArray.size(); i++) {
        if (removeSet.count(hdr->samples[i]) == 1) {
            removeArray[i] = true;
            n_removed++;
        }
    }
    std::cerr << n_removed << " individuals removed" << std::endl;
    return n_removed;
}

int readDiscardSet(
    const char *discard_fin,
    std::set<std::string> &discardSet
) {
    if (!discard_fin) return 0;

    std::ifstream discardFile(discard_fin);
    assert(discardFile);
    std::string var;

    while (discardFile >> var) {
        discardSet.insert(var);
    }
    std::cerr << "Pruning " << discardSet.size() << " variants" << std::endl;
    return discardSet.size();
}

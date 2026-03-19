#include "hmm.h"
#include "prepare_bcf.h"
#include "genome_rules.h"
#include <boost/program_options.hpp>
#include <htslib/synced_bcf_reader.h>
#include <algorithm>
#include <sstream>

#define MHC_START 27518932
#define MHC_END 33480487
#define KIR_START 54071493
#define KIR_END 54992731

namespace po=boost::program_options;

struct hmm_params_t {
    double phred_minor;
    double phred_major;
    double phred_phase_switch;
    double phred_telomere;
    double phred_err;
    bool unpenalized_phase_switch_between_ps_boundaries;
    std::vector<double> inv_allelic_bias;
    int merge_max_gap_hets;
    double merge_max_vaf_frac;
};

// Default inverse allelic bias values used when none are provided on the CLI
static const std::vector<double> DEFAULT_INV_ALLELIC_BIAS = {
    4.0, 5.0, 6.0, 8.0, 10.0, 15.0, 20.0, 30.0, 50.0,
    80.0, 120.0, 145.0, 170.0, 230.0, 310.0, 400.0, 500.0, 610.0
};
static const std::string DEFAULT_INV_ALLELIC_BIAS_STR = "4,5,6,8,10,15,20,30,50,80,120,145,170,230,310,400,500,610";

static std::string join_values(const std::vector<double> &values) {
    std::ostringstream oss;
    for (size_t i = 0; i < values.size(); ++i) {
        if (i > 0) oss << ' ';
        oss << values[i];
    }
    return oss.str();
}

static std::vector<double> parse_inv_allelic_bias(const std::string &input) {
    std::string normalized = input;
    std::replace(normalized.begin(), normalized.end(), ',', ' ');
    std::stringstream ss(normalized);
    std::vector<double> values;
    double v;
    while (ss >> v) values.push_back(v);
    if (values.empty()) values = DEFAULT_INV_ALLELIC_BIAS; // fall back to default if nothing parsed
    return values;
}

/**
 * @brief compute sample sd from vector
 * 
 * @param v vector of data
 * @param mean mean of v
 * @return double sample sd of v
 */
static double compute_sample_sd(const std::vector<double> &v, double mean) {
    double s2 = 0.0;
    int n = v.size();
    for (int i = 0; i< n; i++) {
        s2 += (v[i] - mean) * (v[i] - mean);
    }
    s2 /= (n-1);
    return std::sqrt(s2);
}

/**
 * @brief compute the BAF and BAF SE in a given region
 * 
 * @param obs list of variants with allelic depths and ref bias
 * @param model HMM object with viterbi path 
 * @param start start index (of obs) defining region to compute BAF from
 * @param end end index (of obs) defining region to compute BAF until
 * @return deviation_t the BAF dev and BAF SE between start and end 
 */
static deviation_t get_deviation(
	const std::vector<variant_t> &obs, 
	const HMM &model, 
	int start,
    int end
) {
    bool phase_switch;
    double bias_factor;
    double tot_denom = 0.0;
    double tot_num = 0.0;

    double denom;
    double num;
    std::vector<double> v_theta;
    for (int i=start; i<=end; i++) {
        // skip variants not in an event
        if (model.path[i] == model.base_state) continue;

        // account for phase switches
        phase_switch = model.path[i] < model.base_state;

        // adjust the observed ADs by the ref bias
        bias_factor = (1-obs[i].bias)/obs[i].bias;
        denom = obs[i].refAD*bias_factor+obs[i].altAD;

        // match phase of variant with phase of the viterbi path
        if (obs[i].phase != phase_switch) {
            num = bias_factor*obs[i].refAD;
        } else {
            num = obs[i].altAD;
        }

        // only use variants with more than 10 fragments 
        if (denom > 10) v_theta.push_back(num/denom);

        tot_num+=num;
        tot_denom += denom;
    }
    int n = v_theta.size();

    double theta = tot_num / tot_denom;

    // (i'm not sure this is the best way to calculate SE but the ref bias adjustment makes things weird)
    double se = compute_sample_sd(v_theta, theta) / std::sqrt(n);

    deviation_t res = {theta, se};
 
    return res;
}

/**
 * @brief extract calls from the viterbi path of the HMM for a given individual
 * 
 * @param obs list of variants for computing BAF deviation
 * @param calls list of calls for individual in question
 * @param model hmm object that has already been fit to the data
 */
static void compute_segment_lods(
    const std::vector<variant_t> &obs,
    std::vector<event_t> &calls,
    HMM &model
) {
    int counter = 0;
    bool in_event;

    int start;
    int end;
    deviation_t deviation;
    double alt;
    double base;
    double lod;

    int n_obs = obs.size();

    for (int i=0; i < n_obs; i++) {
        // logic is a bit confusing but basically we want to get the start and 
        // end indexes of a given segment
        in_event = model.path[i] != model.base_state;
        if (in_event) {
            counter ++;
            if (i != n_obs - 1) {
                continue;
            }
        }

        if (counter == 0) {
            continue;
        }
        
        // one before the start of the call
        start = (in_event) ? i - counter : i - counter - 1;
        // one after the end of the call
        end = (in_event) ? i + 1 : i;

        // get observed deviation and SE in non-zero segment of the call
        deviation = get_deviation(obs, model, start+1, end-1);
        int n_flips = 0;
        alt = model.eval_log_prob(deviation.theta, obs, start, end, &n_flips);
        int n_flips_base = 0;
        base = model.eval_log_prob(0.5, obs, start, end, &n_flips_base);
        lod = (alt - base) / std::log(10);

        event_t event = {(uint32_t) start+1, (uint32_t) end-1, lod, deviation, n_flips};

        calls.push_back(event);
        counter = 0;
    }
}

/**
 * @brief merge adjacent calls with similar VAFs and small gaps
 * 
 * @param obs list of variants with allelic depths
 * @param calls list of calls to merge
 * @param model HMM object with viterbi path
 * @param max_gap_hets maximum number of het sites allowed between calls to merge
 * @param max_vaf_frac maximum fractional VAF difference allowed between calls to merge
 */
static void merge_adjacent_calls(
    const std::vector<variant_t> &obs,
    std::vector<event_t> &calls,
    HMM &model,
    int max_gap_hets = 10,
    double max_vaf_frac = 0.10
) {
    if (calls.size() <= 1) return;
    
    std::vector<event_t> merged_calls;
    merged_calls.push_back(calls[0]);
    
    for (size_t i = 1; i < calls.size(); i++) {
        event_t &prev_call = merged_calls.back();
        event_t &curr_call = calls[i];
        
        // calculate gap between calls in terms of het sites
        int gap_hets = curr_call.start - prev_call.end - 1;
        
        // calculate VAF difference (absolute deviation from 0.5)
        double prev_vaf = std::abs(prev_call.deviation.theta - 0.5);
        double curr_vaf = std::abs(curr_call.deviation.theta - 0.5);
        double vaf_diff = std::abs(prev_vaf - curr_vaf);
        double vaf_scale = std::max(prev_vaf, curr_vaf);
        bool vaf_close = (vaf_scale == 0.0) ? true : (vaf_diff <= vaf_scale * max_vaf_frac);

        // merge if gap is small and VAFs are similar (fractionally)
        if (gap_hets <= max_gap_hets && vaf_close) {
            // extend the previous call to include the current call
            prev_call.end = curr_call.end;
            
            // recalculate deviation and LOD for the merged call
            prev_call.deviation = get_deviation(obs, model, prev_call.start, prev_call.end);
            int n_flips_dummy = 0;
            double alt = model.eval_log_prob(prev_call.deviation.theta, obs, prev_call.start - 1, prev_call.end + 1, &n_flips_dummy);
            int n_flips_base = 0;
            double base = model.eval_log_prob(0.5, obs, prev_call.start - 1, prev_call.end + 1, &n_flips_base);
            prev_call.lod = (alt - base) / std::log(10);

            // update n_flips based on viterbi path and add extra if phase switch between calls
            int n_flips = prev_call.n_flips + curr_call.n_flips;
            if (model.path[prev_call.start] == model.base_state || model.path[curr_call.start] == model.base_state) {
                // keep n_flips as-is if either call starts in base state (shouldn't happen for events)
            } else if (model.path[prev_call.end] == model.get_phase_switch_state(model.path[curr_call.start])) {
                n_flips += 1;
            }
            prev_call.n_flips = n_flips;
        } else {
            // cannot merge, add current call as a new call
            merged_calls.push_back(curr_call);
        }
    }
    
    calls = merged_calls;
}

void write_header(ostream &os) {
    os << "ID" << "\t"; 
    os << "sex" << "\t"; 
    os << "chr" << "\t"; 
    os << "bpStart" << "\t"; 
    os << "bpEnd" << "\t"; 
    os << "length" << "\t"; 
    os << "p" << "\t"; 
    os << "q" << "\t"; 
    os << "nHets" << "\t"; 
    os << "avgHetDist" << "\t"; 
    os << "bdev" << "\t"; 
    os << "bdevSE" << "\t"; 
    os << "depth" << "\t"; 
    os << "depthSE" << "\t"; 
    os << "LOD" << "\t"; 
    os << "nFlips" << endl; 
}

/**
 * @brief write calls for all individuals and a given chromosome
 * 
 * @param os output stream
 * @param calls_map vector containing calls for each individual
 * @param obs list of het variants for each individual (used to report stats)
 * @param hdr input bcf header (used for stats)
 * @param rid 0-indexed chromosome (assuming ordered contigs)
 * @param genome_rules rules annotating chromosome centromeres and short arms
 */
void write_calls(
    ostream &os, 
    std::vector<std::vector<event_t>> &calls_map, 
    std::vector<std::vector<variant_t>> &obs,
    bcf_hdr_t *hdr,
    int rid,
    genome_rules_t *genome_rules
) {
    char p_arm;
    char q_arm;

    // trying to match MoChA's output as much as possible (a few columns have been cut)

    for (int i=0; i < bcf_hdr_nsamples(hdr); i++) {
        int num_calls = calls_map[i].size();
        for (int j=0; j<num_calls; j++) {
            int call_start = obs[i][calls_map[i][j].start].pos;
            int call_end = obs[i][calls_map[i][j].end].pos;
            int call_length = call_end - call_start + 1;
            int n_hets = calls_map[i][j].end - calls_map[i][j].start + 1;
            double avg_het_dist = (double) call_length / (double) n_hets;

            if (calls_map[i][j].start == 0 && !genome_rules->is_short_arm[rid]) {
                p_arm = 'T';
            } else if (call_start < genome_rules->cen_beg[rid]) {
                p_arm = 'Y';
            } else {
                p_arm = 'N';
            }
            if (calls_map[i][j].end == obs[i].size() - 1) {
                q_arm = 'T';
            } else if (call_end > genome_rules->cen_end[rid]) {
                q_arm = 'Y';
            } else {
                q_arm = 'N';
            }

            os << hdr->samples[i] << "\t"; // sample ID
            os << "U" << "\t"; // sex (to be annotated in post)
            os << "chr" << rid+1 << "\t"; // chrom
            os << call_start << "\t"; // start 
            os << call_end << "\t"; // end 
            os << call_length << "\t"; // length
            os << p_arm << "\t"; // p_arm
            os << q_arm << "\t"; // p_arm
            os << n_hets << "\t"; // n_hets
            os << avg_het_dist << "\t"; // het average distance
            os << calls_map[i][j].deviation.theta - 0.5 << "\t"; // deviation from 0.5
            os << calls_map[i][j].deviation.se << "\t"; // s.d. of deviation from 0.5
            os << "nan" << "\t"; // relative depth (to be annotated in post)
            os << "nan" << "\t"; // relative depth sd (to be annotated in post)
            os << calls_map[i][j].lod << "\t"; // lod 
            os << calls_map[i][j].n_flips << endl; // n_flips
        }
        obs[i].clear();
        calls_map[i].clear();
    }
}

/**
 * @brief save the allelic shift from the viterbi path for a given individual
 * 
 * @param model model that has already been fit to data
 * @param sample_obs vector of variants with phase information
 * @param sample_allelic_shift map of position to allelic shift
 */
void save_allelic_shift(
    const HMM &model, 
    const std::vector<variant_t> &sample_obs, 
    std::map<uint32_t, bool> &sample_allelic_shift
) {
    int as;
    int n_sample_obs = sample_obs.size();
    for (int i = 0; i < n_sample_obs; i++) {
        // skip variants that are not in an event
        if (model.path[i] == model.base_state) continue;

        // 1 if shift towards reference, -1 if shift towards alt (i think but it might be flipped)
        as = (model.path[i] > model.base_state) ? 1 : -1;
        as *= (sample_obs[i].phase) ? 1 : -1;
        sample_allelic_shift[sample_obs[i].pos] = (as > 0) ? true : false;
    }
}

/**
 * @brief make calls on a given chromosome for a given individual
 * 
 * @param rid 0-indexed chromosome (assuming ordered contigs)
 * @param sample_obs list of a given individuals het variants
 * @param calls list of calls to update
 * @param sample_allelic_shift map to store the allelic shift 
 * @param hmm_params hmm params
 * @param genome_rules rules annotating chromosome centromeres and short arms
 * @return 0 if successful, 1 if error
 */
int make_calls(
    int rid,
    const std::vector<variant_t> &sample_obs, 
    std::vector<event_t> &calls, 
    std::map<uint32_t, bool> &sample_allelic_shift,
    hmm_params_t hmm_params,
    genome_rules_t *genome_rules
) {
    double phred_minor = hmm_params.phred_minor;
    double phred_major = hmm_params.phred_major;
    double phred_phase_switch = hmm_params.phred_phase_switch;
    double phred_telomere = hmm_params.phred_telomere;
    double phred_err = hmm_params.phred_err;
    const std::vector<double> &inv_allelic_bias = hmm_params.inv_allelic_bias.empty()
        ? DEFAULT_INV_ALLELIC_BIAS
        : hmm_params.inv_allelic_bias;

    // the resulting calls will be stored internally in the model object
    HMM model(
        rid,
        phred_minor,
        phred_major,
        phred_phase_switch,
        phred_telomere,
        phred_err,
        hmm_params.unpenalized_phase_switch_between_ps_boundaries,
        inv_allelic_bias,
        genome_rules
    );
    if (model.viterbi(sample_obs) < 0) {
        return 1;
    }
    
    // each individual has it's own list of calls so we can use threads safely
    compute_segment_lods(sample_obs, calls, model);
    
    // merge adjacent calls with similar VAFs and small gaps
    merge_adjacent_calls(sample_obs, calls, model, hmm_params.merge_max_gap_hets, hmm_params.merge_max_vaf_frac);
    
    save_allelic_shift(model, sample_obs, sample_allelic_shift);
    return 0;
}

/**
 * @brief write the output bcf header
 * 
 * @param bcf_out output bcf 
 * @param hdr input bcf header
 * @return bcf_hdr_t* output bcf header with a new AS field
 */
static bcf_hdr_t *print_hdr(htsFile *bcf_out, bcf_hdr_t *hdr) {
    bcf_hdr_t *hdr_out = bcf_hdr_dup(hdr);
    int as_id = bcf_hdr_id2int(hdr, BCF_DT_ID, "AS");
    if(!bcf_hdr_idinfo_exists(hdr, BCF_HL_FMT, as_id)) {
        bcf_hdr_append(
            hdr_out,
            "##FORMAT=<ID=AS,Number=1,Type=Integer,Description=\"Allelic shift (1/-1 if the alternate allele is over/under represented)\">"
        );
    }
    if(bcf_hdr_write(bcf_out, hdr_out) < 0) {
        cerr << "Unable to write to output VCF" << endl;
        exit(1);
    }
    return hdr_out;
}

/**
 * @brief function to write allelic shift direction at each het variant 
 * 
 * @param allelic_shift_map vector containing map from het variant to allelic shift for each individual
 * @param rid 0-indexed chromosome (assuming contigs are in order)
 * @param sr input bcf
 * @param hdr input bcf header 
 * @param bcf_out output bcf
 * @param hdr_out output bcf header
 * @return int 0 if successful, 1 if unsuccessful
 */
int write_allelic_shift_bcf(
    std::vector<std::map<uint32_t, bool>> &allelic_shift_map, 
    int rid,
    bcf_srs_t *sr,
    bcf_hdr_t *hdr,
    htsFile *bcf_out,
    bcf_hdr_t *hdr_out
) {
    bcf_sr_seek(sr, bcf_hdr_id2name(hdr, rid), 0);
    int nsmpl = bcf_hdr_nsamples(hdr_out);
    int *allelic_shift = (int *) calloc(nsmpl, sizeof(int));
    while ( bcf_sr_next_line(sr) ) {    
        bcf1_t *rec = bcf_sr_get_line(sr, 0);
        uint32_t pos = rec->pos+1;
        if(rec->rid != rid) break;
        for (int i = 0; i < nsmpl; i++) {
            // 0 if site is not within a call
            if (allelic_shift_map[i].count(pos) == 0) {
                allelic_shift[i] = 0;
                continue;
            }
            // 1 if towards reference, -1 if away from reference
            allelic_shift[i] = (allelic_shift_map[i][pos]) ? 1 : -1;
        }
        bcf_update_format_int32(hdr_out, rec, "AS", allelic_shift, (int) nsmpl);
        if (bcf_write(bcf_out, hdr_out, rec) < 0) {
            cerr << "Unable to write to output BCF" << endl;
            return 1;
        }
    }
    free(allelic_shift);
    return 0;
}

/**
 * @brief read the allelic depths at het sites skipping over masked regions
 * 
 * @param sr input bcf
 * @param hdr input bcf header
 * @param rid 0-indexed chromosome to read from
 * @param start start bp to read from
 * @param end end bp to read until
 * @param obs a list of het variants for each sample
 * @param refBiasMap ref bias at each variant
 * @param cnvLookup map with CNVs called from depth
 * @param genome_rules rules annotating chromosome centromeres and short arms
 * @return int the number of variants included in the model
 */
int read_bcf_ad(
    bcf_srs_t *sr, 
    bcf_hdr_t *hdr, 
    int rid, 
    uint32_t start,
    uint32_t end,
    vector<vector<variant_t>> &obs, 
    map<string, double> &refBiasMap,
    std::map<string, std::queue<CNV>> &cnvLookup,
    genome_rules_t *genome_rules
) {
    int *GT = NULL, GT_buf_size = 0;
    int *AD = NULL, AD_buf_size = 0;
    int *PS = NULL, PS_buf_size = 0;

    char variant_name[1000]; 

    uint8_t refAD;
    uint8_t altAD;
    bool phase;
    int32_t phase_set;
    double bias;
    int n = 0;

    // seek to the correct chromosome and starting location
    bcf_sr_seek(sr, bcf_hdr_id2name(hdr, rid), start);
    int nsmpl = bcf_hdr_nsamples(hdr);

    // obs will have a vector of variants for each individual
    obs.resize(nsmpl);

    std::vector<std::queue<CNV>> cnvLookupArray(nsmpl);
    for (int i = 0; i < nsmpl; i++) {
        cnvLookupArray[i] = cnvLookup[hdr->samples[i]];
    }

    while ( bcf_sr_next_line(sr) ) {
        bcf1_t *rec = bcf_sr_get_line(sr, 0);
        if (rec->rid != rid) break;
        if (rec->pos > end) break;

        // skip variants that have been blacklisted 
        if(variantInCNV(rid+1, rec->pos+1, cnvLookup["BLACKLIST"])) continue;
        assert(rec->n_allele == 2);
        assert(bcf_unpack(rec, BCF_UN_STR) == 0);

        assert(bcf_get_format_int32(hdr, rec, "AD", &AD, &AD_buf_size) == 2*nsmpl);
        assert(bcf_get_genotypes(hdr, rec, &GT, &GT_buf_size) == 2*nsmpl);
        int ps_ret = bcf_get_format_int32(hdr, rec, "PS", &PS, &PS_buf_size);
        int ps_vals_per_sample = (ps_ret > 0 && ps_ret % nsmpl == 0) ? ps_ret / nsmpl : 0;

        // skip variants in MHC and KIR
        if ((rid+1 == 6) && (rec->pos+1 > MHC_START) && (rec->pos+1 < MHC_END)) continue; 
        if ((rid+1 == 19) && (rec->pos+1 > KIR_START) && (rec->pos+1 < KIR_END)) continue; 

        // apply individual level CNV masks on the read ADs
        maskCNV(rid+1, rec->pos+1, hdr, cnvLookupArray, AD);

        snprintf(variant_name, sizeof(variant_name), "%d_%ld_%s_%s", rid+1, rec->pos+1, rec->d.allele[0], rec->d.allele[1]);

        // if the variant is missing from refBias, we give it a value of -1 (will be skipped)
        bias = (refBiasMap.count(variant_name) == 1) ? refBiasMap[variant_name] : -1;
        if (bias > 0.55 || bias < 0.45) continue;
        assert(bias > 0);

        // skip variants in short arms
        if(genome_rules->is_short_arm[rid] && rec->pos+1 < genome_rules->cen_beg[rid]) continue;

        // skip variants in the centromeres
        if(rec->pos+1 > genome_rules->cen_beg[rid] && rec->pos+1 < genome_rules->cen_end[rid]) continue;

        n++;
        for (int i = 0; i < nsmpl; i++) {
            //form the vector of het variants
            if (AD[2*i] == bcf_int32_missing || AD[2*i+1] == bcf_int32_missing) continue;
            if (bcf_gt_is_missing(GT[2*i]) || bcf_gt_is_missing(GT[2*i+1])) continue;
            if (!bcf_gt_is_phased(GT[2*i+1])) continue;
            refAD = (uint8_t) AD[2*i];
            altAD = (uint8_t) AD[2*i+1];
            phase = (bool) bcf_gt_allele(GT[2*i+1]);
            phase_set = PHASE_SET_MISSING;
            if (ps_vals_per_sample > 0 && PS != NULL) {
                int ps_idx = i * ps_vals_per_sample;
                if (PS[ps_idx] != bcf_int32_missing && PS[ps_idx] != bcf_int32_vector_end) {
                    phase_set = PS[ps_idx];
                }
            }

            variant_t var = {refAD, altAD, phase, phase_set, (uint32_t) rec->pos+1, bias};
            obs[i].push_back(var);
        }
    }
    
    free(GT);
    free(AD);
    free(PS);

    return n;
}

/**
 * @brief make calls for all individuals on a single chromosome
 * 
 * @param obs a list of het variants for each sample
 * @param rid 0-indexed chromosome number to process
 * @param sr synced reader object for input bcf 
 * @param hdr bcf header
 * @param call_os output stream to write calls to
 * @param bcf_out output bcf to write allelic shift to
 * @param hdr_out output bcf header
 * @param hmm_params parameters for the hmm
 * @param genome_rules rules annotating chromosome centromeres and short arms
 */
void process_chrom( 
    std::vector<std::vector<variant_t>> &obs, 
    int rid, 
    bcf_srs_t *sr,
    bcf_hdr_t *hdr, 
    ostream &call_os, 
    htsFile *bcf_out,
    bcf_hdr_t *hdr_out,
    hmm_params_t hmm_params,
    genome_rules_t *genome_rules
) {
    int nsmpl = bcf_hdr_nsamples(hdr);
    std::vector<std::vector<event_t>> calls_map;
    std::vector<std::map<uint32_t, bool>> allelic_shift_map;
    calls_map.resize(nsmpl);
    allelic_shift_map.resize(nsmpl);

    cerr << "Done reading chromosome " << rid+1 << endl;
    if (genome_rules->is_short_arm[rid]) {
        cerr << "Ignoring variants in the short arm" << endl;
    }
    #pragma omp parallel for
    for (int i=0; i < nsmpl; i++) {
        std::vector<event_t> sample_calls;
        std::map<uint32_t, bool> sample_allelic_shift;

        if(make_calls(rid, obs[i], sample_calls, sample_allelic_shift, hmm_params, genome_rules)) {
            cerr << "Warning: no variants for " << hdr->samples[i] << endl;
        }

        allelic_shift_map[i] = sample_allelic_shift;
        calls_map[i] = sample_calls;
    }
    write_calls(call_os, calls_map, obs, hdr, rid, genome_rules);
    write_allelic_shift_bcf(allelic_shift_map, rid, sr, hdr, bcf_out, hdr_out);
    cerr << "Done making calls for chromosome " << rid+1 << endl;
}

int make_calls_main(std::vector<std::string> &argv) {
    const char *ref_bias_fin = NULL;
    const char *cnv_mask_fin = NULL;
    const char *calls_fout = NULL;
    const char *bcf_fout = NULL;
    const char *bcf_fin = NULL;
    
    int rChr0=-1, rStart=0, rEnd=1e9;
    bool mask_deletions=false;


    po::variables_map vm;
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "produce help message")
        ("mask-deletions", "flag to set whether or not to mask deletions called from depth")
        ("ref-bias", po::value<std::string>(), "file to read ref-bias from")
        ("cnv-mask", po::value<std::string>(), "file to read cnv-mask from")
        ("calls-out", po::value<std::string>(), "file to write calls to")
        ("bcf", po::value<std::string>(), "bcf to read from")
        ("bcf-out", po::value<std::string>(), "bcf file to write HMM path to")
        ("region", po::value<std::string>(), "region to run HMM on")
        ("phred-minor", po::value<double>()->default_value(40.0), "HMM minor transition probability (phred)")
        ("phred-major", po::value<double>()->default_value(50.0), "HMM major transition probability (phred)")
        ("phred-phase-switch", po::value<double>()->default_value(40.0), "HMM phase switch probability (phred)")
        ("unpenalized-phase-switch-ps-boundary", "Do not penalize phase-switch transitions that occur at phase-set boundaries (PS tag changes)")
        ("phred-telomere", po::value<double>()->default_value(20.0), "HMM starting in dev state penalty (phred)")
        ("phred-error", po::value<double>()->default_value(15.0), "HMM crop threshold (phred)")
        ("inv-allelic-bias", po::value<std::string>()
            ->default_value(DEFAULT_INV_ALLELIC_BIAS_STR, "4,5,6,8,10,15,20,30,50,80,120,145,170,230,310,400,500,610"),
            "Comma-separated inverse allelic bias values (default mirrors existing tuning)")
        ("merge-max-gap-hets", po::value<int>()->default_value(10), "Maximum number of het sites between calls to merge")
        ("merge-max-vaf-frac", po::value<double>()->default_value(0.10), "Maximum fractional VAF difference between calls to merge");

	try {

		po::store(po::command_line_parser(argv).options(desc).run(), vm);
		po::notify(vm);

    } catch(exception& e) {
        cerr << "Error: " << e.what() << endl;
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
    if (vm.count("bcf")) {
        cerr << "Reading bcf from " << vm["bcf"].as<std::string>() << endl;
        bcf_fin = vm["bcf"].as<std::string>().c_str();
    } else {
        cerr << "Error: need to supply a bcf" << endl;
        return 1;
    }	

    if (vm.count("ref-bias")) {
        cerr << "Reading ref bias from " << vm["ref-bias"].as<std::string>() << endl;
        ref_bias_fin = vm["ref-bias"].as<std::string>().c_str();
    } else {
        cerr << "No ref bias supplied" << endl;
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
        cerr << "No CNV mask supplied" << endl;
    }	

    if (vm.count("calls-out")) {
        cerr << "Writing calls to " << vm["calls-out"].as<std::string>() << endl;
        calls_fout = vm["calls-out"].as<std::string>().c_str();
    } else {
        cerr << "Writing calls to stdout" << endl;
    }	

    if (vm.count("bcf-out")) {
        cerr << "Writing allelic shift to " << vm["bcf-out"].as<std::string>() << endl;
        bcf_fout = vm["bcf-out"].as<std::string>().c_str();
    } else {
        cerr << "Writing allelic shift to stdout" << endl;
        bcf_fout = "-";
    }	

    if (vm.count("region")) {
        string regionStr = vm["region"].as<std::string>();
        cerr << "Running HMM on region: " << regionStr << endl;
        int colonPos = regionStr.find(':');
        assert(colonPos != (int) string::npos);
        rChr0 = parseChrStrToChr0(regionStr.substr(0, colonPos).c_str());
        assert(sscanf(regionStr.substr(colonPos+1).c_str(), "%d-%d", &rStart, &rEnd) == 2);
        rStart --;
        rEnd--;
        assert(rStart >= 0);
        assert(rEnd > rStart);
    } else {
        cerr << "Running HMM on entire bcf file" << endl;
        rChr0 = -1;
        rStart = 0;
        rEnd = 1e9;
    }

    hmm_params_t hmm_params = {
        vm["phred-minor"].as<double>(),
        vm["phred-major"].as<double>(),
        vm["phred-phase-switch"].as<double>(),
        vm["phred-telomere"].as<double>(),
        vm["phred-error"].as<double>(),
        vm.count("unpenalized-phase-switch-ps-boundary") > 0,
        parse_inv_allelic_bias(vm["inv-allelic-bias"].as<std::string>()),
        vm["merge-max-gap-hets"].as<int>(),
        vm["merge-max-vaf-frac"].as<double>()
    };

    cerr << "HMM minor phred: " << hmm_params.phred_minor << "\n";
    cerr << "HMM major phred: " << hmm_params.phred_major << "\n";
    cerr << "HMM phase switch phred: " << hmm_params.phred_phase_switch << "\n";
    cerr << "Unpenalized phase switches at PS boundaries: " << (hmm_params.unpenalized_phase_switch_between_ps_boundaries ? "true" : "false") << "\n";
    cerr << "HMM telomere phred: " << hmm_params.phred_telomere << "\n";
    cerr << "HMM error phred: " << hmm_params.phred_err << "\n";
    cerr << "Inverse allelic bias values: " << join_values(hmm_params.inv_allelic_bias) << "\n";
    cerr << "Merge max gap hets: " << hmm_params.merge_max_gap_hets << "\n";
    cerr << "Merge max VAF frac: " << hmm_params.merge_max_vaf_frac << "\n";


    bcf_srs_t *sr = bcf_sr_init();
    bcf_sr_set_opt(sr, BCF_SR_REQUIRE_IDX);
    bcf_sr_set_opt(sr, BCF_SR_PAIR_LOGIC, BCF_SR_PAIR_EXACT);

    if(!bcf_sr_add_reader(sr, bcf_fin)) {
        cerr << "Error: failed to read bcf file " << bcf_fin << endl;
        return 1;
    }
    bcf_hdr_t *hdr = bcf_sr_get_header(sr, 0); assert(hdr != NULL);

    genome_rules_t *genome_rules;
    char build[] = "GRCh38";
    genome_rules = genome_init_alias(stderr, build, hdr); 
    readlist_short_arms(genome_rules, "13,14,15,21,22,chr13,chr14,chr15,chr21,chr22", hdr);

    // loading refBias to calibrate the emission probabilities
    map<string, double> refBiasMap;
    if (ref_bias_fin) {
        ifstream refBiasIn;
        refBiasIn.open(ref_bias_fin); assert(refBiasIn);
        refBiasMap = parseRefBias(&refBiasIn);
        refBiasIn.close();
    }

    // loading CNVs to skip over when reading in het sites
    map<string, queue<CNV>> cnvLookup;
    if (cnv_mask_fin) {
        ifstream germlineCNVsFile;
        germlineCNVsFile.open(cnv_mask_fin); assert(germlineCNVsFile);
        cnvLookup = parseCNVs(&germlineCNVsFile, mask_deletions);
        germlineCNVsFile.close();

        for (int i = 0; i < bcf_hdr_nsamples(hdr); i++) {
            cerr << "Number of CNVs masked in ";
            cerr << hdr->samples[i] << ": ";
            cerr << cnvLookup[hdr->samples[i]].size() << endl;
        }
        cerr << "Number of regions blacklisted: " << cnvLookup["BLACKLIST"].size() << endl;
    }

    std::ostream *call_os = &cout;
    std::ofstream call_file;

    if (calls_fout) {
        call_file.open(calls_fout); assert(call_file);
        call_os = &call_file;
    } 

    htsFile *bcf_out = hts_open(bcf_fout, "wb");
    bcf_hdr_t *hdr_out = print_hdr(bcf_out, hdr);

    std::vector<std::vector<variant_t>> obs;

    // this if statment might be unnecessary but it allows us to run the HMM on 
    // a specific focal region in case we want to be more generous with the 
    // prior
    if (rChr0 > 0 ) {
        assert(read_bcf_ad(sr, hdr, rChr0, rStart, rEnd, obs, refBiasMap, cnvLookup, genome_rules) > 0);
        write_header(*call_os);
        process_chrom(obs, rChr0, sr, hdr, *call_os, bcf_out, hdr_out, hmm_params, genome_rules);
    } else {
        // process each chromosome in the input file separately
        // most input bcfs for our pipeline contain just 1 chromosome so this 
        // loop should end up skipping the chromosomes not present 
        write_header(*call_os);
        for (int rid = 0; rid < hdr->n[BCF_DT_CTG]; rid++) {
            if(read_bcf_ad(sr, hdr, rid, rStart, rEnd, obs, refBiasMap, cnvLookup, genome_rules) == 0) continue;
            process_chrom(obs, rid, sr, hdr, *call_os, bcf_out, hdr_out, hmm_params, genome_rules);
        }
    }

    if (calls_fout) call_file.close();
    bcf_sr_destroy(sr);
    bcf_hdr_destroy(hdr_out);
    hts_close(bcf_out);
    genome_destroy(genome_rules);

    return 0;
}

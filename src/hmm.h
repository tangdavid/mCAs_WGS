#ifndef __HMM_H__
#define __HMM_H__

#include <vector>
#include <stdint.h>
#include <cmath>
#include <string>
#include <limits>
#include "genome_rules.h"

static const int32_t PHASE_SET_MISSING = std::numeric_limits<int32_t>::min();

struct deviation_t {
    double theta;
    double se;
};

struct variant_t {
    uint8_t refAD;
    uint8_t altAD;
    bool phase;
    int32_t phase_set;
    uint32_t pos;
    double bias;
};

struct event_t {
    uint32_t start;
    uint32_t end;
    double lod;    
    deviation_t deviation; 
    int n_flips;
};


class HMM {
public:
    std::vector<int> path;
    int rid;
    double log_p_minor;
    double log_p_major;
    double log_p_phase_switch;
    double log_p_telomere;
    double log_p_base;
    double log_p_base_extend;
    double log_p_event_extend;
    double err_log_p;
    int n_states;
    int base_state;
    bool unpenalized_phase_switch_at_ps_boundary;
    std::vector<double> allelic_bias;
    genome_rules_t *genome_rules;

    HMM(
        int chrom_id,
        double phred_minor, 
        double phred_major, 
        double phred_phase_switch,
        double phred_telomere,
        double phred_err,
        bool unpenalized_phase_switch_between_ps_boundaries,
        const std::vector<double> &inv_allelic_bias,
        genome_rules_t *gr
    ) :
        rid(chrom_id),
        log_p_minor(-phred_minor/10 * std::log(10)),
        log_p_major(-phred_major/10 * std::log(10)),
        log_p_phase_switch(-phred_phase_switch/10 * std::log(10)),
        log_p_telomere(-phred_telomere/10 * std::log(10)),
        err_log_p(-phred_err/10 * std::log(10)),
        unpenalized_phase_switch_at_ps_boundary(unpenalized_phase_switch_between_ps_boundaries),
        genome_rules(gr)
    {
        /*
        log_p_base = std::log(1 - std::exp(log_p_telomere) * (n_states - 1));
        log_p_base_extend = std::log(1 - std::exp(log_p_event_start) * (n_states - 1));
        log_p_event_extend = std::log(1 - std::exp(log_p_phase_switch) - std::exp(log_p_event_stop));
        */

        n_states = inv_allelic_bias.size() * 2 + 1; 
        log_p_base = 0;
        log_p_base_extend = 0;
        log_p_event_extend = 0;

        base_state = n_states/2;
        
        allelic_bias.resize(n_states);
        allelic_bias[base_state] = 0.5;
        for (int s = 0; s < n_states/2; s++) {
            allelic_bias[base_state - s - 1] = 0.5 - 1/inv_allelic_bias[s];
            allelic_bias[base_state + s + 1] = 0.5 + 1/inv_allelic_bias[s];
        }
    }
    double get_log_transition(int s0, int s1, bool in_p_arm, bool phase_set_boundary=false);
    int get_phase_switch_state(int state);
    bool is_phase_set_boundary(const std::vector<variant_t> &obs, int idx);
    double log_emit(int state, variant_t var);
    double log_emit_fixed_theta(double theta, variant_t var);
    int viterbi(const std::vector<variant_t> &obs); 
    double eval_log_prob(double theta, const std::vector<variant_t> &obs,int start,int end, int *n_flips);
};

#endif

#include "hmm.h"
#include <cfloat>

/**
 * @brief rescale log probabilities to prevent overflow 
 * (subtract the max log probability from each entry)
 * 
 * @param log_prb 
 */
static inline void rescale_log_prb(std::vector<double> &log_prb) {
    double max = -DBL_MAX;
    int n = log_prb.size();
    for (int i = 0; i < n; i++) max = (max > log_prb[i]) ? max : log_prb[i];
    for (int i = 0; i < n; i++) log_prb[i] -= max;
}

/**
 * @brief compute the viterbi path given a set of variant ADs with phase and ref bias
 * 
 * @param obs list of variants with phase and ref bias
 * @return int -1 if no observations, 0 if successful
 */
int HMM::viterbi(const std::vector<variant_t> &obs) {

    int n_obs = obs.size();

    if (n_obs == 0) return -1;
    std::vector<std::vector<double>> dp(n_obs, std::vector<double> (n_states));
    std::vector<std::vector<int>> backtrack(n_obs, std::vector<int> (n_states));
    path.resize(n_obs, 0);
   
    for (int s0 = 0; s0 < n_states; s0++) {
        dp[0][s0] = HMM::log_emit(s0, obs[0]);
        if (s0 != base_state) {
            // should we give telomere boost to short arm chromosomes? 
            // the telomere boost is applied to the first and last observation 
            // regardless of whether or not they are actually in the telomere
            // dp[0][s0] += (genome_rules->is_short_arm[rid]) ? log_p_minor : log_p_telomere;
            dp[0][s0] += log_p_telomere;
        } else {
            dp[0][s0] += log_p_base;
        }
    }

    for (int i = 1; i < n_obs; i++) {
        bool phase_set_boundary = HMM::is_phase_set_boundary(obs, i);
        for (int s1 = 0; s1 < n_states; s1++) {
            int max_idx = 0;  
            double max_val = -DBL_MAX;
            double current_val;
            for (int s0 = 0; s0<n_states; s0++) {
                // disallow changes between non-zero states
                if (
                    s0 != base_state && 
                    s1 != base_state && 
                    s0 != s1 && 
                    s0 != HMM::get_phase_switch_state(s1)
                ) {
                    continue;
                }
                current_val = dp[i-1][s0] + HMM::get_log_transition(s0, s1, obs[i].pos < genome_rules->cen_beg[rid], phase_set_boundary);
                if (current_val > max_val) {
                    max_val = current_val;
                    max_idx = s0;
                }
            }
            dp[i][s1] = max_val + HMM::log_emit(s1, obs[i]);
            backtrack[i][s1] = max_idx;
        }
        // normalize probabilities every 1000 iterations to prevent overflow
        if (i % 10000 ) rescale_log_prb(dp[i]);
    }

    int max_idx = 0;
    double max_val = -DBL_MAX;
    // add closing costs to HMM and store the max index for path decoding
    for (int s0 = 0; s0 < n_states; s0++) {
        if (s0 != base_state) {
            dp[n_obs-1][s0] += log_p_telomere;
        } else {
            dp[n_obs-1][s0] += log_p_base;
        }
        if (dp[n_obs-1][s0] > max_val) {
            max_val = dp[n_obs-1][s0];
            max_idx = s0;
        }
    }
    
    // backtrack to get the viterbi path
    path[n_obs-1] = max_idx;
    for (int i = n_obs-1; i >0; i--) {
        path[i-1] = backtrack[i][path[i]];
    }

    return 0;
}

/**
 * @brief get transition probabilities from state s0 to state s1
 * 
 * @param s0 
 * @param s1 
 * @param in_p_arm 
 * @return double 
 */
double HMM::get_log_transition(int s0, int s1, bool in_p_arm, bool phase_set_boundary) {

    // flip the start and stop probabilities in q arm for symmetry
    double log_p_event_start = in_p_arm ? log_p_major : log_p_minor;
    double log_p_event_stop = in_p_arm ? log_p_minor : log_p_major;

    if (s0 != base_state && s1 == base_state) {
        return log_p_event_stop;
    }
    else if (s0 == base_state && s1 != base_state) {
        return log_p_event_start;
    }
    else if (s0 == s1 && s0 != base_state) {
        return log_p_event_extend;
    }
    else if (s0 == s1 && s0 == base_state) {
        return log_p_base_extend;
    }
    else if (s0 == HMM::get_phase_switch_state(s1)) {
        if (unpenalized_phase_switch_at_ps_boundary && phase_set_boundary) {
            return 0;
        }
        return log_p_phase_switch;
    }
    else {
        // this handles impossible transitions 
        // maybe we could have small boosts from transitioning between non-zero
        // states
        return -DBL_MAX;
    }
}

bool HMM::is_phase_set_boundary(const std::vector<variant_t> &obs, int idx) {
    if (idx <= 0 || idx >= (int) obs.size()) return false;
    int32_t prev_phase_set = obs[idx-1].phase_set;
    int32_t curr_phase_set = obs[idx].phase_set;
    if (prev_phase_set == PHASE_SET_MISSING || curr_phase_set == PHASE_SET_MISSING) return false;
    return prev_phase_set != curr_phase_set;
}

/**
 * @brief get the phase switch state for a given state
 * states are coded such that the n_states/2 is the base state and deviations 
 * above and below correspond to positive vs negative shifts towards the ref alt
 *
 * @param state 
 * @return int 
 */
int HMM::get_phase_switch_state(int state) {
    return n_states - state - 1;
}


static inline double log_choose(int n, int k) {
    return std::lgamma(double(n+1)) - std::lgamma(double(k+1)) - std::lgamma(double(n-k+1));
}

static inline double log_binomial_pmf(double p, int n, int k) {
    return log_choose(n, k) + double(k)*std::log(p) + double(n-k)*std::log(1-p);
}

/**
 * @brief emission probabilities from a given state (using internal theta values)
 * 
 * @param state 
 * @param var 
 * @return double 
 */
double HMM::log_emit(int state, variant_t var) {
    return HMM::log_emit_fixed_theta(allelic_bias[state], var);
}

/**
 * @brief emission probabilities from a given state (with a passed theta)
 * 
 * @param state 
 * @param var 
 * @param theta the fraction of reads from hap1
 * @return double 
 */
double HMM::log_emit_fixed_theta(double theta, variant_t var) {
    // phase == true when ref is on hap1
    // phase == false when ref is on hap2
    // theta > 0.5 is when hap1 is represented more than hap2
    // theta < 0.5 is when hap2 is represented more than hap1
    // i.e. theta is fraction of reads from hap1

    // number of observed ref reads is drawn from binomial(n, p)
    double p;
    if ( var.phase ) { 
        // ref count is from hap1 
        p=(theta/(1-theta))/(theta/(1-theta) + ((1-var.bias)/var.bias));
    } else {
        // ref count is from hap2
        p=((1-theta)/theta)/((1-theta)/theta + ((1-var.bias)/var.bias));
    }

    //double min_emit_thresh = log_binomial_pmf(p, 50, std::round(50 * p)) + err_log_p;
    double min_emit_thresh = -2.1 + err_log_p;
    double emit_prob = log_binomial_pmf(p, var.refAD + var.altAD, var.refAD);

    // crop outliers
    return (emit_prob < min_emit_thresh) ? min_emit_thresh : emit_prob;
}

/**
 * @brief compute log probability given a fixed value of theta (for LOD calculations)
 * 
 * @param theta fraction of reads from hap1
 * @param obs variants with phase, ADs, and ref bias
 * @param start index of observation 1 before the start of the segment of interest
 * @param end index of observation 1 after the end of the segment of interest
 * @param n_flips to keep track of number of flips in the segment 
 * @return double 
 */
double HMM::eval_log_prob(
    double theta,
    const std::vector<variant_t> &obs,
    int start,
    int end,
    int *n_flips
) {
    // assumes that theta is greater than 0.5
    theta = (theta > 0.5) ? theta : 1 - theta;
    double res = 0;
    *n_flips = 0;
    int n_penalized_flips = 0;
    double t;
    variant_t var;
    int n_obs = obs.size();

    if (start == -1) {
        //res += (path[start+1] == base_state) ? log_p_base : log_p_telomere;
        start += 1;
    }

    if (end == n_obs) {
        //res += (path[end-1] == base_state) ? log_p_base : log_p_telomere;
        end -= 1;
    }

    for (int i = start+1; i <= end; i++) {
        // looks like MoChA's lod calculations do not include the event start 
        // and end penalties
        // uncomment out the block below to add an event start penalty in lod 
        // calculations

        /* 
        res += (theta != 0.5) ? 
            HMM::get_log_transition(path[i-1], path[i])
            : HMM::get_log_transition(base_state, base_state);
        */
        if (path[i-1] == HMM::get_phase_switch_state(path[i])) {
            *n_flips += 1;
            if (!(unpenalized_phase_switch_at_ps_boundary && HMM::is_phase_set_boundary(obs, i))) {
                n_penalized_flips += 1;
            }
        }   
    }
    res += n_penalized_flips * log_p_phase_switch;

    for (int i = start; i <= end; i++) {
        var = obs[i];

        if (path[i] > base_state) {
            t = theta;
        } else if (path[i] < base_state) {
            t = 1 - theta;
        } else {
            t = 0.5;
        }

        res += HMM::log_emit_fixed_theta(t, var);
    }

    return res;
}

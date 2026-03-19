import numpy as np
import subprocess
import gzip
import argparse
import sys
import scipy.stats

def get_allele_frequencies(file_path, num_lines=1000, seed=None):
    np.random.seed(seed)
    AF = []
    with gzip.open(file_path, 'rt', encoding='utf-8') as file:
        for i, line in enumerate(file):  # Start index at 1
            freq = float(line.split('\t')[4])
            freq = np.minimum(freq, 1-freq)
            if freq < 1e-5: continue

            if len(AF) < num_lines:
                AF.append(freq)
            else:
                idx = np.random.randint(0, i)
                if idx < num_lines:
                    AF[idx] = freq
            
    AF = np.array(AF)
    return AF

def get_allele_frequencies_fast(file_path, num_lines=1000):
    """Uses shuf via subprocess to quickly shuffle and read a gzipped file."""
    cmd = f"zcat {file_path} | tail +2 | shuf -n {num_lines}" 
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True, check=True)
    AF = np.array([float(line.split('\t')[4]) for line in result.stdout.splitlines()])
    AF = np.minimum(AF, 1-AF)
    AF = np.maximum(AF, 1/4e5)
    return AF

def simulate_trait(AF, alpha, h2 = 1, n=10000, prevalence=0.002):
    haps = np.random.binomial(1, AF, size = (n, 2, AF.shape[0]))
    beta = np.random.normal(loc = 0, scale = np.sqrt(2* (AF*(1-AF))**alpha))
    G = (haps - haps.mean(axis = (0, 1))) @ beta
    epsilon = np.random.normal(loc = 0, scale = np.sqrt((1-h2)*np.var(G[:,0] - G[:, 1])/h2), size = len(G))
    prolif_score = (G[:, 0 ] - G[:, 1] + epsilon)

    thresh = np.quantile(np.abs(prolif_score), 1-prevalence)
    allelic_shift = (prolif_score>thresh).astype(int) - (prolif_score<-thresh).astype(int)
    return haps, beta, G, prolif_score, allelic_shift

def simulate_GWAS(haps, allelic_shift):
    hets = (haps.sum(axis = 1) == 1)
    cases = allelic_shift != 0 
    a = (cases[:, None] & hets).sum(axis = 0)
    b = ((~cases[:, None]) & hets).sum(axis = 0)
    c = ((cases[:, None]) & (~hets)).sum(axis = 0)
    d = ((~cases[:, None]) & (~hets)).sum(axis = 0)
    fisher_p = (1 - scipy.stats.hypergeom.cdf(a-1, a+b+c+d, a+b, a+c))

    het_on_hap1 = haps[:,0,:] > haps[:,1,:]
    het_on_hap2 = haps[:,1,:] > haps[:,0,:]
    var_hom = het_on_hap1[allelic_shift>0].sum(axis = 0) + het_on_hap2[allelic_shift<0].sum(axis = 0)
    var_rem = het_on_hap2[allelic_shift>0].sum(axis = 0) + het_on_hap1[allelic_shift<0].sum(axis = 0)
    binom_p = scipy.stats.binom(p=0.5,n=var_hom+var_rem).cdf(np.minimum(var_hom, var_rem)) * 2

    return fisher_p, binom_p

def evaluate_power(AF, fisher_p, binom_p, sig_thresh=1e-4):
    bins = np.array([1e-5, 1e-4, 1e-3, 1e-2, 1])
    counts, _ = np.histogram(AF, bins = bins)
    fisher_hits, _ = np.histogram(AF[np.where(fisher_p<sig_thresh)[0]], bins = bins)
    binom_hits, _ = np.histogram(AF[np.where(binom_p<sig_thresh)[0]], bins = bins)

    power = np.zeros((len(bins)-1, 3))
    power[:,0] = counts
    power[:,1] = fisher_hits/np.maximum(1, counts)
    power[:,2] = binom_hits/np.maximum(1, counts)

    return bins, power


def run_simulation(alpha, h2, num_causal, prevalence, repeats=30, AF_path=None, n=10000):
    res = []
    for i in range(repeats):
        print(f"Repeat: {i+1}/{repeats}", file=sys.stderr, flush=True)
        if AF_path is None:
            AF = get_allele_frequencies_fast('/mnt/project/lohdata/david/resources/chr21.afreq.txt.gz', num_lines=num_causal)
        else:
            AF = get_allele_frequencies(AF_path, num_lines=num_causal, seed=None)
        haps, _, _, _, allelic_shift  = simulate_trait(AF, alpha=alpha, h2=h2, n=n)
        fisher_p, binom_p = simulate_GWAS(haps, allelic_shift)
        bins, power = evaluate_power(AF, fisher_p, binom_p)
        res.append(power)

    power = np.array(res)
    mean_power = np.mean(power, axis = 0)
    stderr_power = np.std(power, axis = 0)/np.sqrt(repeats)
    print("ALPHA\tH2\tNUM_CAUSAL\tPREVALENCE\tAF\tNUM_VARIANTS\tFISHER_POWER\tFISHER_POWER_STDERR\tBINOM_POWER\tBINOM_POWER_STDERR")
    for af_bin, mean, stderr in zip(bins,mean_power, stderr_power):
        print(
            alpha,
            h2,
            num_causal,
            prevalence,
            af_bin, 
            np.round(mean[0], 2), 
            np.round(mean[1], 6), 
            np.round(stderr[1], 6), 
            np.round(mean[2], 6), 
            np.round(stderr[2], 6), 
            sep='\t'
        )

if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("alpha", type=float)
    parser.add_argument("h2", type=float)
    parser.add_argument("num_causal", type=int)
    parser.add_argument("prevalence", type=float)
    parser.add_argument("--n", type=int, default=50000)
    parser.add_argument("--AF", type=str)
    args = parser.parse_args()
    print(f"""Running simulation with: 
    alpha={args.alpha}
    h2={args.h2}
    number_of_causal_var={args.num_causal}
    prevalence={args.prevalence}
    n={args.n}
    """, file=sys.stderr)
    run_simulation(args.alpha, args.h2, args.num_causal, args.prevalence, n=args.n, AF_path=args.AF, repeats=50)

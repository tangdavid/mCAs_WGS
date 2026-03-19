import pandas as pd
import numpy as np
import scipy.stats
import sys

def main():
    df = pd.read_csv(sys.stdin, sep='\t', header=None)
    ref_exome = df.iloc[:, 2]
    alt_exome = df.iloc[:, 3]
    z_exome = df.iloc[:, 4]
    wgs_bdev = df.iloc[:, 5]
    wgs_bdev_se = df.iloc[:, 6]
    n = ref_exome + alt_exome

    expected_p = wgs_bdev + 0.5
    expected_mu = expected_p * n
    expected_sd = np.sqrt(n*expected_p*(1 - expected_p) + n*(n-1)*wgs_bdev_se**2)
    
    expected_replication = [
        1 - scipy.stats.norm.cdf(0.5 * n[i], loc = expected_mu[i], scale=expected_sd[i])
        if n[i] > 0
        else
        0
        for i in range(len(n))
    ]
    expected_replication = np.array(expected_replication, dtype=np.float32)

    df['expected_replication'] = expected_replication
    df.to_csv(sys.stdout, sep='\t', index=None, header=None)
    
if __name__=='__main__':
    main()


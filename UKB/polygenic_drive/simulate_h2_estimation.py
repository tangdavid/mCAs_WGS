import numpy as np
import argparse
import struct

def simulate(N=1000, M=2000, h2=0.3, seed=None):
    np.random.seed(seed)
    p = np.random.beta(1, 1, M)
    X = np.random.binomial(p=p, n=1, size=(N, 2, M)).astype(np.float64)
    for i in range(N):
        thresh = int(np.random.uniform(0, N))
        X[:,0,:thresh]  = 0
        X[:,1,:thresh] = 0
        #print(X)
    X = X[:, 0, :] - X[:, 1, :]
    X /= np.sqrt(2*p*(1-p))
    beta = np.random.normal(size=M)
    epsilon = np.random.normal(size=N)
    G = X @ beta
    G /= np.std(G)
    G = np.sqrt(h2) * G
    epsilon = np.sqrt((1-h2)) * epsilon
    Y = G + epsilon
    K = 1/M * X @ X.T
    K /= np.diag(K).mean()
    Y = (Y > 0).astype(np.uint8)
    return Y, X, K 


def write_gcta_input(K, Y, output_prefix):
    # Write grm.bin (lower triangular part of K)
    N = K.shape[0]
    with open(f"{output_prefix}.grm.bin", "wb") as f:
        for i in range(N):
            for j in range(i + 1):
                f.write(struct.pack("f", K[i, j]))

    # Write grm.id (individual IDs)
    with open(f"{output_prefix}.grm.id", "w") as f:
        for i in range(N):
            f.write(f"{i+1}\t{i+1}\n")

    # Write pheno file
    with open(f"{output_prefix}.pheno", "w") as f:
        for i in range(N):
            f.write(f"{i+1}\t{i+1}\t{Y[i]}\n")

if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Simulate data for h2 estimation")
    parser.add_argument("--N", type=int, default=200, help="Number of individuals")
    parser.add_argument("--M", type=int, default=2000, help="Number of SNPs")
    parser.add_argument("--h2", type=float, default=0.3, help="Heritability")
    parser.add_argument("--seed", type=int, default=None, help="Random seed")
    parser.add_argument("--out-dir", type=str, default=".", help = "Output dir")
    args = parser.parse_args()

    Y, X, K = simulate(M=args.M, N=args.N, h2=args.h2, seed=args.seed)
    write_gcta_input(K, Y, f"{args.out_dir}/{args.N}_{args.M}_{args.h2}_{args.seed}")

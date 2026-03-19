import numpy as np
import pandas as pd
import sys
import scipy.optimize
import argparse

def compute_overlap(df_ID):
    n = len(df_ID)
    overlap = np.zeros((n, n))
    i = 0
    for _,r1 in df_ID.iterrows():
        j = 0
        for _,r2 in df_ID.iterrows():
            start = max(r1['bpStart'], r2['bpStart'])
            end = min(r1['bpEnd'], r2['bpEnd'])
            overlap[i, j] = (end - start)/1e3
            if r1['bpStart'] >= r2['bpEnd']: overlap[i, j] = 0
            if r1['bpEnd'] <= r2['bpStart']: overlap[i, j] = 0
            j+=1
        i +=1
    return overlap

def var_selection_one_round(overlap, n, readsMissing, expectedReads, mask, orientation):
    best_error = np.inf
    lb = np.array([0 if x=="OUTER" else -np.inf for x in orientation])
    ub = np.array([0 if x=="INNER" else np.inf for x in orientation])

    for i in range(n):
        if mask[i]: continue
        mask[i] = True
        # depthDev = np.linalg.solve(overlap[:, mask].T @ overlap[:, mask], overlap[:, mask].T @ readsMissing)
        #depthDev = -depthDev.reshape(-1,1)

        solution = scipy.optimize.lsq_linear(overlap[:, mask], readsMissing.reshape(-1,), bounds=(lb[mask], ub[mask]))
        depthDev = solution.x.reshape(-1,1)

        squared_error = np.sum((overlap[:, mask] @ depthDev - readsMissing)**2/expectedReads)

        if squared_error < best_error:
            best_region = i
            best_depth_dev = depthDev
            best_error = squared_error
        mask[i] = False
    mask[best_region] = True
    rescale = expectedReads.max()/overlap.max()
    #print(rescale, expectedReads[0]/overlap[0][0], file=sys.stderr, flush=True)
    return best_error, best_depth_dev/rescale

def var_selection(df_ID, zscore_cutoff):
    n = len(df_ID)
    overlap = compute_overlap(df_ID)
    expectedReads = df_ID['EXPinclCNV'].to_numpy().reshape(-1, 1)
    readsMissing = (df_ID['OBSinclCNV'] - df_ID['EXPinclCNV']).to_numpy().reshape(-1, 1)
    squared_error = float('inf')
    
    orientation = df_ID['orientation'].to_numpy()
    mask = np.array([True if x=="SINK" else False for x in orientation])

    solution = scipy.optimize.lsq_linear(overlap[:, mask], readsMissing.reshape(-1,))
    depthDev = solution.x.reshape(-1,1)
    squared_error = np.sum((overlap[:, mask] @ depthDev - readsMissing)**2/expectedReads)

    num_regions=mask.sum()
    while(squared_error > zscore_cutoff**2 and mask.sum() < n):
        squared_error, depthDev = var_selection_one_round(overlap, n, readsMissing, expectedReads, mask, orientation)
        #print(n, mask, file=sys.stderr, flush=True)
        assert mask.sum() == num_regions+1
        num_regions=mask.sum()
        # print(squared_error)
        # print(mask)
        # print(depthDev)
    res = np.hstack([
        df_ID.iloc[mask, [0,1,2,3]].to_numpy(), 
        depthDev.round(4),
        np.repeat(squared_error.round(4), len(depthDev)).reshape(-1, 1),
    ])
    nonzeros = (depthDev.round(4) != 0).flatten()
    res = res[nonzeros]
    return res


def resolve_overlapping(df, zscore_cutoff):
    df = df.dropna()
    df['zscore'] = (df['OBSinclCNV']/df['EXPinclCNV'] - 1)/(np.sqrt(df['EXPinclCNV'])/df['EXPinclCNV'])
    sink = (df['orientation']=="SINK") 
    inner = (df['orientation']=="INNER") 
    outer = (df['orientation']=="OUTER") 
    df = df[sink | inner | outer]
    groups = df.groupby(['ID', 'chr'])
    res = []
    for ID, df_ID in groups:
        mat = var_selection(df_ID, zscore_cutoff)
        res.append(mat)
    res = pd.DataFrame(np.vstack(res))
    res.columns = ['ID', 'chr', 'bpStart', 'bpEnd', 'depthDev', 'unexplainedResidual']
    res['depthDev'] = res['depthDev'].astype(float)
    res = res.merge(df, on=['ID', 'chr', 'bpStart', 'bpEnd'])
    res = res.infer_objects()
    return res
    


if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--zscore', default=3, type=float)
    args = parser.parse_args()
    df = pd.read_csv(sys.stdin, sep='\t')
    res = resolve_overlapping(df=df, zscore_cutoff=args.zscore)
    res.to_csv(sys.stdout, index=False, sep ='\t', float_format='%.4f')

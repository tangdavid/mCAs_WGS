import numpy as np
import pandas as pd
import re
import sys

def assign_event_type(df):
    X = df['bdev'].to_numpy()
    Y = df['depth'].to_numpy()
    assignment = np.repeat('CN-LOH', len(df))

    vdj_regions = [
        'chr2:85000000-92000000',
        'chr7:35000000-40000000',
        'chr7:140000000-145000000',
        'chr14:20000000-25000000',
        'chr14:102000000-108000000',
        'chr22:20000000-25000000'
    ]
    vdj_query = ' or '.join(
        [
            f'(chr=="{chrom}" and bpStart>{bpStart} and bpEnd<{bpEnd})' 
            for chrom,bpStart,bpEnd in [
                re.split(':|-', r) 
                for r in vdj_regions
            ]
        ]
    )
    vdj_events_index = df.reset_index().query(vdj_query).index.to_numpy()
    assignment[vdj_events_index] = 'LOSS'
    gain_separator = 0.8*X
    loss_separator = -0.6*X
    assignment[Y>gain_separator] = 'GAIN'
    assignment[Y<loss_separator] = 'LOSS'

    return assignment

def filter_germline_events(df):
    num_events = len(df)
    bsize = 1e6
    df['bin'] = [
        f'{chrom}:{start:0.0f}-{end:0.0f}' 
        for chrom,start,end 
        in zip(
            df['chr'], 
            round(df['bpStart']/bsize)*bsize, 
            round(df['bpEnd']/bsize)*bsize
        )
    ]
    x = df \
        .query('type=="GAIN"') \
        .groupby('bin').agg(
            count=('type', len),
            cf_mean=('cf', np.mean),
            cf_std=('cf', np.std),
            len_mean=('length', np.mean),
            len_std=('length', np.std)
        ) \
        .query('count > 5 and len_mean < 1e7 and cf_mean > 0.3 and cf_std < cf_mean * 0.3') \
        .reset_index() 
    germline_filter = set(x['bin'])
    df = df.query('not(bin in @germline_filter and type == "GAIN")')
    df = df.query('not(chr=="chr8" and bpStart>86e6 and bpEnd<86.5e6)') # looks germline from defining focal regions
    print(f"Filtered {num_events - len(df)} putatively germline duplications", file=sys.stderr)
    print("Regions:", germline_filter, file=sys.stderr)
    num_events=len(df)
    df = df.query('not(chr=="chr4" and bpStart<40e6 and bpEnd>39e6 and type!="GAIN" and length<2e6)') # repeat expansion at RFC1
    print(f"Filtered {num_events - len(df)} events and RFC1", file=sys.stderr)
    num_events=len(df)
    df = df.query('not(chr=="chr10" and bpStart<112e6 and bpEnd>111e6 and type!="GAIN" and length<2e6)') # repeat expansion at FRA10B
    print(f"Filtered {num_events - len(df)} events at FRA10B", file=sys.stderr)
    num_events=len(df)
    df = df.drop('bin', axis = 1)
    return df

def filter_tri21(df):
    num_events = len(df)
    df_chrom_depth = pd.read_csv('/mnt/project/lohdata/resources/WGS_depth/PCadj/per_chr_depth/ID_40709.per_chr_depth_PCadj.txt.gz', sep='\t')
    trisomy21_IDs = set(df_chrom_depth.query('chr21>1.4')['ID'])
    res = df[[(not (ID in trisomy21_IDs and chrom=="chr21")) for ID,chrom in zip(df['ID'], df['chr'])]].reset_index(drop=True)
    print(f"Filtered {num_events - len(res)} mCAs from {len(trisomy21_IDs)} individuals with tri21", file=sys.stderr)
    return res

def filter_sample_contamination(df):
    num_events = len(df)
    df_sample_contamination = pd.read_csv('/mnt/project/lohdata/david/mCAs_WGS/sample_data/ID.sample_contamination.txt', sep='\t')
    blacklist = set(df_sample_contamination.query('sample_contamination>0.02')['ID'])
    res = df[[(not ID in blacklist) for ID in df['ID']]].reset_index(drop=True)
    print(f"Filtered {num_events - len(res)} events from {len(blacklist)} individuals with sample contamination > 0.02", file=sys.stderr)
    return res

def filter_interstitial_CNLOH(df):
    num_events = len(df)
    interstitial_CNLOH_per_ind = df\
    .query('type=="CN-LOH" and p!="T" and q!="T"')\
    .groupby('ID')\
    .agg(
        count = ('ID', len)
    )\
    .reset_index()
    blacklist = set(interstitial_CNLOH_per_ind.query('count>50')['ID'])
    res = df[[(not ID in blacklist) for ID in df['ID']]].reset_index(drop=True)
    print(f"Filtered {num_events - len(res)} events from {len(blacklist)} individuals with >50 interstitial CN-LOH calls", file=sys.stderr)
    return res

if __name__=="__main__":
    df = pd.read_csv("/mnt/project/lohdata/david/mCAs_WGS/calls/WGS_500k.calls.raw.txt.gz", sep='\t')
    num_events = len(df)
    df_withdrawn=pd.read_csv('/mnt/project/lohdata/resources/BOLT-LMM/w40709_CURRENT.FID_IID.txt', sep=' ', header=None)
    df_withdrawn.columns=['FID', 'IID']
    withdrawn = set(df_withdrawn['FID'])
    df = df.query('not ID in @withdrawn')
    print(f"Filtered {num_events - len(df)} events from withdrawn individuals", file=sys.stderr)

    num_events = len(df)
    df = df.dropna()
    print(f"Filtered {num_events - len(df)} events due to missing depth", file=sys.stderr)

    assignments = assign_event_type(df)
    df['type'] = assignments
    df.loc[df['type'] == 'CN-LOH', 'cf'] = 2 * df['bdev']
    df.loc[df['type'] != 'CN-LOH', 'cf'] = 2 * np.abs(df['depth'])
    df = filter_germline_events(df)
    df = filter_tri21(df)
    df = filter_sample_contamination(df)
    df = filter_interstitial_CNLOH(df)


    df_age = pd.read_csv('/mnt/project/lohdata/david/mCAs_WGS/sample_data/ID_age.40709.txt', sep='\t')
    df_age['ageBin'] = df_age['age']//5 * 5
    df = df.merge(df_age, on = 'ID', how='left')
    df.to_csv(sys.stdout, sep='\t', index=None)

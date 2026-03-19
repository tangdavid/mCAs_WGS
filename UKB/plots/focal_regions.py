import numpy as np
import pandas as pd
from collections import defaultdict


def define_focal_regions(calls):
    df=calls.copy()
    df = df[(df['length']<2e6) & ((df['type']=='GAIN') | (df['type']=='LOSS'))]
    bins = {'LOSS':np.zeros((22, 300),int), 'GAIN':np.zeros((22, 300),int)}
    bins2rows = defaultdict(list)

    blacklist = [
        'chr2:85000000-92000000',
        'chr7:35000000-40000000',
        'chr7:140000000-145000000',
        'chr14:20000000-25000000',
        'chr14:102000000-108000000',
        'chr22:20000000-25000000',
        # 'chr2:23000000-27000000',
        # 'chr4:100000000-110000000',
        # 'chr13:48000000-52000000',
        # 'chr16:29000000-32000000',
        # 'chr17:30000000-33000000',
        # 'chr22:28000000-30000000'
    ]
    blacklist_bins = set()
    for region in blacklist:
        chrom = region.split(':')[0]
        start = int(region.split(':')[1].split('-')[0])/1e6
        end = int(region.split(':')[1].split('-')[1])/1e6
        for i in range(int(start), int(end)+1): blacklist_bins.add((chrom, i))

    for idx, row in df.iterrows():
        for i in range(int(row['bpStart']/1e6), int(row['bpEnd']/1e6)+1):
            if (row['chr'], i) in blacklist_bins: continue
            zero_based_chr_idx = int(row['chr'][3:])-1
            bins[row['type']][zero_based_chr_idx][i] += 1
            bins2rows[(row['type'], zero_based_chr_idx, i)].append(idx)

    keep = list()
    current_region = []
    last_bin = None
    for copy_number in ['LOSS', 'GAIN']:
        for i, chrom_bins in enumerate(bins[copy_number]):
            chrom = i+1
            for j, counts in enumerate(chrom_bins):
                if (len(current_region) == 0): start = j
                if counts>10:
                    current_region += bins2rows[(copy_number, i, j)]
                    continue
                if (len(current_region) == 0): continue
                unique_ids = np.unique(current_region)
                keep += current_region
                df.loc[unique_ids,'region'] = f'chr{chrom}:{start*1e6:0.0f}-{j*1e6:0.0f}'
                current_region.clear()

    keep = np.unique(keep)
    df = df.loc[keep,:]
    return df

def compute_focal_index(df):
    focal_index = {}
    focal_index_position = {}
    for (region, copy_number), df_region in df.groupby(['region', 'type']):
        chrom = region.split(':')[0]
        start = int(region.split(':')[1].split('-')[0])
        end = int(region.split(':')[1].split('-')[1])
        counts = defaultdict(set)
        for _, row in df_region.iterrows():
            for i in range(start, end+1, 1000):
                if i > row['bpStart'] and i < row['bpEnd']:
                    counts[i].add(row['ID'])
        position_value = sorted([(k,len(v)/len(df_region['ID'].unique())) for k,v in counts.items()], key=lambda x:x[1], reverse=True)[0]
        focal_index_position[(region, copy_number)], focal_index[(region, copy_number)] = position_value
    
    return focal_index, focal_index_position

def assign_region_to_locus(df, gtf, cyto, chip_genes, focal_index_position):
    gtf['coding_gene'] = ['NM_' in x for x in gtf['transcript']]
    transcript_table = gtf.query('annot=="transcript"')
    coding_genes = set(transcript_table.query('coding_gene')['gene'])

    gene_density = defaultdict(dict)

    num_events = df.groupby(['region', 'type']).agg(count= ('ID', len)).to_dict()['count']
    for i, row in df.iterrows():
        chrom = row['chr']
        start = row['bpStart']
        end = row['bpEnd']
        key = (row['region'], row['type'])
        for j,transcript in transcript_table.query('chr==@chrom and start<@end and end>@start').iterrows():
            counts = gene_density[key].get(transcript['gene'], 0)
            gene_density[key][transcript['gene']] = counts+1/num_events[key]

    region_to_locus = {}
    for (region,copy_number), position in focal_index_position.items():
        chrom = region.split(':')[0]
        region_to_locus[(region, copy_number)] = f'{chrom[3:]}{cyto[(cyto[0]==chrom) & (cyto[1] <= position) & (cyto[2]>=position)].iloc[0, 3]}'

    for (region, copy_number),gene_counts in gene_density.items():
        if len(chip_genes.intersection(gene_counts.keys())) > 0:
            region_to_locus[(region, copy_number)] = chip_genes.intersection(gene_counts.keys()).pop()
            continue
        gene_counts = sorted([(k, v) for k, v in gene_counts.items() if k in coding_genes], reverse=True, key=lambda x: x[1])
        sorted_gene_counts = [x[1] for x in gene_counts]
        sorted_gene_names = [x[0] for x in gene_counts]
        
        if len(sorted_gene_counts) == 0: continue
        if sorted_gene_counts[0] < 0.5 or (len(sorted_gene_counts)>1 and sorted_gene_counts[1]*2 > sorted_gene_counts[0]): continue
        region_to_locus[(region, copy_number)] = sorted_gene_names[0]
    
    return df.merge(
        pd.DataFrame(
            [(region, cn, locus) for (region, cn), locus in region_to_locus.items()], 
            columns=['region', 'type', 'locus']
        ), 
            on=['region', 'type']
    ), gene_density

def get_snp_array_counts(df, snp_array):
    regions = np.array(
        [
            (chrom[3:], interval.split('-')[0], interval.split('-')[1])
            for chrom,interval in [(x.split(':')) for x in np.unique(df['region'])]
        ]
        , dtype=int)

    for i,row in snp_array.iterrows():
        for region in regions:
            if row['CHR'] != region[0]: continue
            if row['START_MB'] < region[1]/1e6 or row['END_MB'] > region[2]/1e6: continue
            snp_array.loc[i,'region'] = f'chr{region[0]}:{region[1]}-{region[2]}'
    res = snp_array.dropna().groupby('region').agg(
        snp_array_count = ('region', len)
    )
    return res
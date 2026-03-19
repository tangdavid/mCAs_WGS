import numpy as np
import pandas as pd
import sys

def main():
	df_cancer = pd.read_csv(sys.stdin, sep='\t')
	df_ID = pd.read_csv('/mnt/project/lohdata/david/mCAs_WGS/sample_data/WGS_500k.included_IDs.txt')
	df_1Mb = pd.read_csv('/mnt/project/lohdata/david/mCAs_WGS/del13q/WGS_500k.del13q.calls.depth.txt.gz', sep='\t')
	df_hmm = pd.read_csv('/mnt/project/lohdata/david/mCAs_WGS/calls/WGS_500k.calls.txt', sep='\t')

	df_pheno = df_1Mb[['ID', 'del13q_depth', 'depthDev']].merge(
		df_hmm.query('chr=="chr13" and (type=="LOSS" or type=="CN-LOH") and bpEnd>49e6 and bpStart<51e6')\
			.groupby('ID') \
			.agg(
				del13q_baf = ('type', lambda x: np.any(x=='LOSS')),
				cnloh13q = ('type', lambda x: np.any(x=='CN-LOH')),
				cf = ('cf', lambda x: np.max(x))
			),
		on='ID',
		how='outer'
	).merge(df_ID, how='outer').fillna(0).merge(
		df_hmm.query('chr=="chr12" and type=="GAIN" and p=="T" and q=="T"').assign(tri12=1)[['ID', 'tri12', 'cf']],
		how = 'outer',
		on = 'ID'
	).fillna(0).assign(
		del13q = lambda x: ((x['del13q_depth'] == 1) | x['del13q_baf']).astype(int),
		cnloh13q = lambda x: x['cnloh13q'].astype(int),
		tri12 = lambda x: x['tri12'].astype(int),
		cf = lambda x: np.maximum(np.maximum(x['cf_x'], x['cf_y']), -2*x['depthDev']),
		FID=lambda x: x['ID'],
		IID=lambda x: x['ID'],
	).assign(
		del13q_cnloh13q = lambda x: x['del13q'] | x['cnloh13q'],
		del13q_tri12 = lambda x: x['del13q'] | x['tri12'],
		del13q_cnloh13q_tri12 = lambda x: x['del13q'] | x['cnloh13q'] | x['tri12'],
	)[['FID', 'IID', 'cf', 'del13q', 'tri12', 'cnloh13q', 'del13q_cnloh13q', 'del13q_tri12', 'del13q_cnloh13q_tri12']].astype(str)

	df_pheno.merge(df_cancer.astype(str), how="outer").fillna('NA').to_csv(sys.stdout, sep='\t', index=None)

if __name__=="__main__":
	main()

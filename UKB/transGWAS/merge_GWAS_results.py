import numpy as np
import pandas as pd
import gc

def downsample_gwas(df, log10p_col='LOG10P'):
    ps = 10**(np.ceil(df[log10p_col] * 2)/2) * 1e-2
    df = df[(np.random.uniform(size=len(ps)) < ps) | (df[log10p_col]>4)]
    return df

def main():
    GWAS_DIR='/mnt/project/lohdata/david/mCAs_WGS/del13q/GWAS'
    index_variants = pd.read_csv(f'{GWAS_DIR}/C911.CLL.GWAS_catalog.top_pval.txt', sep='\t')[['CHROM', 'POS']]

    for i, pheno in enumerate(['del13q_tri12', 'C911']):
        print(pheno)
        res = []
        for chrom in range(1, 22+1):
            df = pd.read_csv(
                f'{GWAS_DIR}/ukb_step2_CLL_GWAS_c{chrom}_{pheno}.regenie.gz', 
                sep=' ', 
                dtype={'EXTRA': object}
        )
            df = downsample_gwas(df)
            res.append(df)
        df = pd.concat(res).reset_index(drop=True).drop('EXTRA', axis=1)
        del(res)
        gc.collect()

        df['NOVEL'] = False
        thresh = -np.log10(5e-8)
        for var in df.query('LOG10P> @thresh').itertuples():
            df.loc[var.Index, 'NOVEL'] = True
            for index_var in index_variants.itertuples():
                if var.CHROM == index_var.CHROM and abs(var.GENPOS - index_var.POS) < 0.5e6:
                    df.loc[var.Index, 'NOVEL'] = False
                    break
        
        seen = []
        df['INDEX'] = False
        for var in df.query('LOG10P>@thresh').sort_values('LOG10P', ascending=False).itertuples():
            df.loc[var.Index, 'INDEX'] = True
            for index_var in seen:
                if var.CHROM == index_var.CHROM and abs(var.GENPOS - index_var.GENPOS) < 0.5e6:
                    df.loc[var.Index, 'INDEX'] = False
                    break
            seen.append(var)

        df.to_csv(f'{pheno}.GWAS.regenie.downsampled.txt.gz', sep='\t', index=None)

if __name__=="__main__":
    main()

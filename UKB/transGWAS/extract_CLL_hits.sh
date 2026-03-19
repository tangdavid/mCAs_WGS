set -euo pipefail

TMP_DIR=$HOME/tmp
mkdir -p $TMP_DIR

curl -s https://www.ebi.ac.uk/gwas/api/search/downloads/full \
    | tr ' ' '_'  \
    | awk -F'\t' -v OFS='\t' '
NR==1 {
    for(i=1;i<=NF;i++) f[$i]=i; 
} $2==28165464 {
    split($f["STRONGEST_SNP-RISK_ALLELE"], snp, "-"); 
    print $f["CHR_ID"],$f["CHR_POS"],$f["P-VALUE"],$f["PVALUE_MLOG"],snp[1],snp[2],$f["OR_or_BETA"],$f["95%_CI_(TEXT)"]
}' \
    | tr -d '[]' \
    > $TMP_DIR/CLL.GWAS_catalog.txt

ls /mnt/project/lohdata/david/mCAs_WGS/del13q/GWAS/*c22*.regenie.gz \
    | xargs -I{} basename -s .regenie.gz {} \
    | cut -f6- -d_ \
    | while read PHENO; do
    echo "Extracting GWAS catalog hits for $PHENO..." >> /dev/stderr
    awk '
    ARGIND==1 {
        extract[$1":"$2]=$6
        odds_ratio[$1":"$2]=$7;
        or_ci_text[$1":"$2]=$8;
    } 
    ARGIND==2 {
        if (FNR==1) {
            print $0,"CLL_OR","CLL_OR_LOWER","CLL_OR_UPPER";
        }
        if ($1":"$2 in extract) {
            split(or_ci_text[$1":"$2], or_ci, "-");
            if (extract[$1":"$2] == $5) {
                print $0,odds_ratio[$1":"$2], or_ci[1], or_ci[2];
            } else {
                print $0,1/odds_ratio[$1":"$2], 1/or_ci[2], 1/or_ci[1];
            }
        }
    }' $TMP_DIR/CLL.GWAS_catalog.txt <(zcat /mnt/project/lohdata/david/mCAs_WGS/del13q/GWAS/ukb_step2_CLL_GWAS_c{1..22}_${PHENO}.regenie.gz ) \
        | tr ' ' '\t' > $OUT_DIR/$PHENO.CLL.GWAS_catalog.index_var.txt
    awk '
    function abs(x) {return (x>0)?x:-x} 
    ARGIND==1 {locus[$1][$2]=1}
    ARGIND==2 && $1 in locus {
        for(pos in locus[$1]) if(abs($2-pos) < 0.5e6) if(maxp[$1][pos] < $13) maxp[$1][pos]=$13
    } 
    END {
        print "CHROM", "POS", "LOG10P_BEST"
        for(chrom in locus) for(pos in locus[chrom]) print chrom, pos, maxp[chrom][pos]
    }' $TMP_DIR/CLL.GWAS_catalog.txt <(zcat /mnt/project/lohdata/david/mCAs_WGS/del13q/GWAS/ukb_step2_CLL_GWAS_c{1..22}_$PHENO.regenie.gz) \
        | tr ' ' '\t' > $OUT_DIR/$PHENO.CLL.GWAS_catalog.top_pval.txt
done


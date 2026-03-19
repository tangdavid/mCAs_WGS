set -euo pipefail


export TMP_DIR=$HOME/tmp
mkdir -p $TMP_DIR

format_regenie_sumstats () {
    PHENO=$1
    if [ ! -f /mnt/project/lohdata/david/mCAs_WGS/del13q/GWAS/ukb_step2_CLL_GWAS_c1_${PHENO}.regenie.gz ]; then
        echo "No summary stats for $PHENO"
        exit 1
    else
        echo "Munging $PHENO sumstats..."
    fi
    awk '
    BEGIN {print "SNP","CHR","BP","A1","A2","A1FREQ","N","BETA","SE","CHISQ","P"} 
    ARGIND==1 {SNP[$1"-"$2]=$4} 
    ARGIND==2 && FNR==1 {for(i=1;i<=NF;i++) f[$i]=i} 
    ARGIND==2 && $1"-"$2 in SNP {
        print SNP[$1"-"$2],$f["CHROM"],$f["GENPOS"],$f["ALLELE0"],$f["ALLELE1"],1-$f["A1FREQ"],$f["N"],$f["BETA"],$f["SE"],$f["CHISQ"],10^(-$f["LOG10P"])
    }' \
        /mnt/project/lohdata/david/ldsc/weights/weights.hm3_noMHC.bed \
        <(zcat /mnt/project/lohdata/david/mCAs_WGS/del13q/GWAS/ukb_step2_CLL_GWAS_c{1..22}_${PHENO}.regenie.gz) \
        > $TMP_DIR/$PHENO.sumstats
    munge_sumstats.py --sumstats $TMP_DIR/$PHENO.sumstats --out $TMP_DIR/$PHENO.munged
}

format_regenie_sumstats del13q_tri12
format_regenie_sumstats C911

ldsc.py \
    --rg $TMP_DIR/C911.munged.sumstats.gz,$TMP_DIR/del13q_tri12.munged.sumstats.gz \
    --ref-ld-chr /mnt/project/lohdata/david/ldsc/weights/weights.hm3_noMHC. \
    --w-ld-chr /mnt/project/lohdata/david/ldsc/weights/weights.hm3_noMHC. \
    --out $OUT_DIR/C911.del13q_tri12

# set on input:
# OUT_DIR

set -euo pipefail

BOLT_DIR=/mnt/project/lohdata/resources/BOLT-LMM
REGENIE_DIR=/mnt/project/lohdata/david/regenie
NUM_PCS=20

TMP_DIR=$HOME/tmp
mkdir -p $TMP_DIR

pip3 install numpy pandas -q

echo `date`': making phenotype...'
zcat /mnt/project/lohdata/resources/phenotypes/merged_disease_phenos.20250401.tab.gz \
    | awk -v OFS='\t' 'NR==1 {for(i=1;i<=NF;i++) f[$i]=i} {print $f["FID"],$f["IID"],$f["C911"],$f["C911.HES"],$f["C911.cancer"]}'  \
    | python3 make_13q_phenotype.py \
    > $OUT_DIR/CLL_pheno.txt

regenie \
	--step 1 \
	--bed $REGENIE_DIR/ukb_cal_chr1_22_v2_merged \
	--extract $REGENIE_DIR/qc_keep.snplist \
	--remove $REGENIE_DIR/qc_keep.mindrem.id,$BOLT_DIR/bolt.in_plink_but_not_imputed.FID_IID.907.txt,$BOLT_DIR/w40709_CURRENT.FID_IID.txt,$BOLT_DIR/remove.nonEUR.FID_IID.40709.txt \
	--phenoFile $OUT_DIR/CLL_pheno.txt \
	--covarFile $BOLT_DIR/covars.40709.txt.gz \
    --covarCol cov_SEX \
    --covarCol cov_AGE \
    --covarCol cov_AGE_SQ \
    --covarCol cov_SMOKING_STATUS \
    --covarCol PC{1:$NUM_PCS} \
	--bt \
	--bsize 1000 \
	--lowmem \
	--lowmem-prefix $TMP_DIR/regenie_tmp_preds \
    --threads `nproc` \
    --gz \
	--out $OUT_DIR/ukb_step1_mCA_trans_GWAS

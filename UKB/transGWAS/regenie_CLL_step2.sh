# set on input: 
# CHR
# OUT_DIR

set -euo pipefail

TMP_DIR=${HOME}/tmp
mkdir -p $TMP_DIR

BOLT_DIR=/mnt/project/lohdata/resources/BOLT-LMM
REGENIE_DIR=/mnt/project/lohdata/david/regenie
PHENO_DIR=/mnt/project/lohdata/david/mCAs_WGS/del13q/GWAS
BGEN_DIR="/mnt/project/lohdata/resources/TOPMed_common/bgen"
NUM_PCS=20

awk -v PFX=$PHENO_DIR '{
    n=split($2,a,"/"); 
    print $1,PFX"/"a[n]
}' $PHENO_DIR/ukb_step1_mCA_trans_GWAS_pred.list \
    > $TMP_DIR/ukb_step1_mCA_trans_GWAS_pred.list

regenie \
	--step 2 \
	--bgen "$BGEN_DIR/TOPMed_common.chr${CHR}.bgen" \
    --sample "$BGEN_DIR/TOPMed_common.chr${CHR}.sample" \
    --ref-first \
	--remove $REGENIE_DIR/qc_keep.mindrem.id,$BOLT_DIR/bolt.in_plink_but_not_imputed.FID_IID.907.txt,$BOLT_DIR/w40709_CURRENT.FID_IID.txt,$BOLT_DIR/remove.nonEUR.FID_IID.40709.txt \
	--phenoFile $PHENO_DIR/CLL_pheno.txt \
	--covarFile $BOLT_DIR/covars.40709.txt.gz \
    --covarCol cov_SEX \
    --covarCol cov_AGE \
    --covarCol cov_AGE_SQ \
    --covarCol cov_SMOKING_STATUS \
    --covarCol PC{1:$NUM_PCS} \
	--bt \
    --firth --approx --pThresh 0.01 \
	--bsize 400 \
    --pred $TMP_DIR/ukb_step1_mCA_trans_GWAS_pred.list \
    --threads `nproc` \
    --minMAC 1000 \
    --minINFO 0.3 \
    --gz \
	--out $OUT_DIR/ukb_step2_CLL_GWAS_c$CHR

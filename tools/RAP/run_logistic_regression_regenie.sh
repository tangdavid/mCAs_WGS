# set on input: 
# CHR
# OUT_DIR
# PHENO_FILE
# OUT_FILE
# EXCLUDE_FILE

set -euo pipefail

TMP_DIR=${HOME}/tmp
mkdir -p $TMP_DIR

BOLT_DIR=/mnt/project/lohdata/resources/BOLT-LMM
REGENIE_DIR=/mnt/project/lohdata/david/regenie
BGEN_DIR="/mnt/project/lohdata/resources/TOPMed_common/bgen"
NUM_PCS=20

regenie \
	--step 2 \
	--bgen "$BGEN_DIR/TOPMed_common.chr${CHR}.bgen" \
    --sample "$BGEN_DIR/TOPMed_common.chr${CHR}.sample" \
    --ref-first \
	--remove $REGENIE_DIR/qc_keep.mindrem.id,$BOLT_DIR/bolt.in_plink_but_not_imputed.FID_IID.907.txt,$BOLT_DIR/w40709_CURRENT.FID_IID.txt,$BOLT_DIR/remove.nonWhite.FID_IID.40709.txt,$EXCLUDE_FILE \
	--phenoFile $PHENO_FILE \
	--covarFile $BOLT_DIR/covars.40709.txt.gz \
    --covarCol cov_SEX \
    --covarCol cov_AGE \
    --covarCol cov_AGE_SQ \
    --covarCol cov_SMOKING_STATUS \
    --covarCol PC{1:$NUM_PCS} \
	--bsize 400 \
    --bt --firth --approx --pThresh 0.01 \
    --ignore-pred \
    --threads `nproc` \
    --minMAC 1000 \
    --minINFO 0.3 \
    --gz \
	--out $OUT_FILE

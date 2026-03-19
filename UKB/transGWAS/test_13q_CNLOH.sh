# set on input: 
# OUT_DIR

set -euo pipefail

TMP_DIR=${HOME}/tmp
mkdir -p $TMP_DIR

BOLT_DIR=/mnt/project/lohdata/resources/BOLT-LMM
REGENIE_DIR=/mnt/project/lohdata/david/regenie
BGEN_DIR="/mnt/project/Bulk/DRAGEN WGS/DRAGEN population level WGS variants, BGEN format [500k release]/"
NUM_PCS=20

pip3 install numpy pandas -q

echo `date`': making phenotype...'
awk -v OFS='\t' 'BEGIN {print "FID","IID";}' \
    | python3 make_13q_phenotype.py \
    | awk '$4!=1 {print $1,$2,$6}' \
    > $TMP_DIR/13qCNLOH.pheno.txt

regenie \
	--step 2 \
	--bgen "$BGEN_DIR/ukb24309_c13_b0_v1.bgen" \
    --sample "$BGEN_DIR/ukb24309_c13_b0_v1.sample" \
	--remove $REGENIE_DIR/qc_keep.mindrem.id,$BOLT_DIR/bolt.in_plink_but_not_imputed.FID_IID.907.txt,$BOLT_DIR/w40709_CURRENT.FID_IID.txt,$BOLT_DIR/remove.nonEUR.FID_IID.40709.txt \
	--phenoFile $TMP_DIR/13qCNLOH.pheno.txt \
	--covarFile $BOLT_DIR/covars.40709.txt.gz \
    --covarCol cov_SEX \
    --covarCol cov_AGE \
    --covarCol cov_AGE_SQ \
    --covarCol cov_SMOKING_STATUS \
    --covarCol PC{1:$NUM_PCS} \
	--bt \
    --firth --approx --pThresh 0.01 \
	--bsize 400 \
    --ignore-pred \
    --range chr13:48000000-52000000 \
    --extract <(awk '$7=="PASS" && $2>48e6 && $2<52e6 {split($8,a,";"); if(int(substr(a[1],4)) < 30) next; print substr($3,8)}' /mnt/project/Bulk/DRAGEN\ WGS/DRAGEN\ population\ level\ WGS\ variants\,\ PLINK\ format\ \[500k\ release\]/ukb24308_c13_b0_v1.pvar ) \
    --threads `nproc` \
    --minMAC 30 \
    --minINFO 0.3 \
    --gz \
	--out $OUT_DIR/13qSNPs

# set on input: 
# CHR
# OUT_DIR

set -euo pipefail

TMP_DIR=${HOME}/tmp
mkdir -p $TMP_DIR

BOLT_DIR=/mnt/project/lohdata/resources/BOLT-LMM
REGENIE_DIR=/mnt/project/lohdata/david/regenie
PHENO_DIR=/mnt/project/lohdata/david/mCAs_WGS/del13q/CLL_progression
BGEN_DIR="/mnt/project/lohdata/resources/TOPMed_common/bgen"
NUM_PCS=20

#zcat $BOLT_DIR/covars.40709.txt.gz | awk -v OFS='\t' 'ARGIND==1{cf[$1]=$3} ARGIND==2 {print $0, (FNR==1) ? "cf": cf[$1]+0}' $PHENO_DIR/CLL_progression_pheno.txt - > $TMP_DIR/covar.txt 

wget 'https://github.com/rgcgithub/regenie/releases/download/v4.1/regenie_v4.1.gz_x86_64_ubuntu20_mkl.zip' -O regenie_v4.1.gz_x86_64_ubuntu20_mkl.zip
unzip regenie_v4.1.gz_x86_64_ubuntu20_mkl.zip

./regenie_v4.1.gz_x86_64_ubuntu20_mkl \
	--step 2 \
	--bgen "$BGEN_DIR/TOPMed_common.chr${CHR}.bgen" \
    --sample "$BGEN_DIR/TOPMed_common.chr${CHR}.sample" \
    --ref-first \
	--remove $REGENIE_DIR/qc_keep.mindrem.id,$BOLT_DIR/bolt.in_plink_but_not_imputed.FID_IID.907.txt,$BOLT_DIR/w40709_CURRENT.FID_IID.txt,$BOLT_DIR/remove.nonEUR.FID_IID.40709.txt \
	--phenoFile <(awk 'ARGIND==1 {pairs[$1][$2]=1; pairs[$2][$1]=1} ARGIND==2 {if($1 in pairs) for(ID in pairs[$1]) if (ID in keep) next; print $0; keep[$1]=1}' "/mnt/project/Bulk/Genotype Results/Genotype calls/ukb_rel.dat" /mnt/project/lohdata/david/mCAs_WGS/del13q/CLL_progression/CLL_progression_pheno.txt) \
	--covarFile <(zcat $BOLT_DIR/covars.40709.txt.gz | awk -v OFS='\t' 'ARGIND==1{cf[$1]=$3} ARGIND==2 {print $0, (FNR==1) ? "cf": cf[$1]+0}' $PHENO_DIR/CLL_progression_pheno.txt -) \
    --covarCol cov_SEX \
    --covarCol cov_AGE \
    --covarCol cov_AGE_SQ \
    --covarCol cov_SMOKING_STATUS \
    --covarCol PC{1:$NUM_PCS} \
    --covarCol cf \
	--bt \
    --firth --approx --pThresh 0.1 \
    --phenoColList cancer \
	--bsize 400 \
    --ignore-pred \
    --threads `nproc` \
    --minMAC 5 \
    --minINFO 0.3 \
    --gz \
	--out $OUT_DIR/chr$CHR.bt.survival

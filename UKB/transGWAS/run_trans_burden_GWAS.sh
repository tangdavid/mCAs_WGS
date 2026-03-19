# set on input: 
# CHR
# OUT_DIR

set -euo pipefail

TMP_DIR=${HOME}/tmp
mkdir -p $TMP_DIR

BOLT_DIR=/mnt/project/lohdata/resources/BOLT-LMM
REGENIE_DIR=/mnt/project/lohdata/david/regenie
BGEN_DIR=/mnt/project/lohdata/nolan/virome/burden_mask_bgen
NUM_PCS=20

cat /mnt/project/lohdata/david/mCAs_WGS/calls/WGS_500k.calls.txt  \
    | awk -v OFS='\t' -v CHR=chr$CHR '
ARGIND==1 && $2!="WGS7" && $2!="WGS8" && $2!="WGS9" {
    samples[$1]=1
}
ARGIND==2 && FNR==1 { for(i=1;i<=NF;i++) f[$i]=i }
ARGIND==2 && $f["chr"] != CHR {
    ID=$(f["ID"])
    type=$(f["type"])
    mca[ID][type]=1
    mca[ID]["ANY"]=1
}
ARGIND==2 {any_mca[$(f["ID"])]=1}
END {
    n=split("FID,IID,CN-LOH,GAIN,LOSS,ANY", cols, ",")
    for (i=1; i<=n; i++) printf "%s%s", cols[i], (i<n? OFS : ORS)
    for (sample in samples) {
        for (i=1; i<=n; i++) printf "%s%s", (i<=2 ? sample : mca[sample][cols[i]]+0), (i<n? OFS : ORS)
    }
}
' /mnt/project/lohdata/resources/WGS_depth/sample_batches/ID_40709.WGS_group.batch.txt - \
    > $TMP_DIR/chr$CHR.mCA_pheno.txt 

plink2 \
    --bgen $BGEN_DIR/chr${CHR}.LoF.scaffold.phased.AS.bgen ref-first \
    --sample $BGEN_DIR/chr${CHR}.LoF.scaffold.phased.AS.fixed.sample \
    --make-pgen erase-phase \
    --out $TMP_DIR/chr${CHR}.LoF

regenie \
	--step 2 \
    --pgen $TMP_DIR/chr${CHR}.LoF \
	--remove $REGENIE_DIR/qc_keep.mindrem.id,$BOLT_DIR/bolt.in_plink_but_not_imputed.FID_IID.907.txt,$BOLT_DIR/w40709_CURRENT.FID_IID.txt,$BOLT_DIR/remove.nonWhite.FID_IID.40709.txt \
	--phenoFile $TMP_DIR/chr$CHR.mCA_pheno.txt \
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
    --threads `nproc` \
    --gz \
    --out $OUT_DIR/chr$CHR.burden_transGWAS

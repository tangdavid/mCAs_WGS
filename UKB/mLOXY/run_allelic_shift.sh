set -euo  pipefail

TMP_DIR=$HOME/tmp
mkdir -p $TMP_DIR

BCF_FILE=/mnt/project/lohdata/david/mCAs_WGS/mLOX_new_phase/WGS_500k.mLOX.allelic_shift.bcf

awk -v CN=$CN '($17==CN && $7=="T" && $8=="T") || NR==1' /mnt/project/lohdata/david/mCAs_WGS/mLOX_new_phase/WGS_500k.mLOX.calls.txt \
    > $TMP_DIR/chrX.$CN.txt 

cat /mnt/project/lohdata/resources/BOLT-LMM/{remove.nonEUR.FID_IID.40709.txt,w40709_CURRENT.FID_IID.txt} \
    | awk '{print $1}' \
    > $TMP_DIR/remove.samples.txt 

BCF_FILE=$BCF_FILE \
CALL_FILE=$TMP_DIR/chrX.$CN.txt \
CHR=chrX \
CN=$CN \
    bash ../cisGWAS/remove_related_individuals.sh >> $TMP_DIR/remove.samples.txt

../../bin/mCAs_WGS compute-cis-fisher \
    --calls $TMP_DIR/chrX.$CN.txt \
    --remove $TMP_DIR/remove.samples.txt \
    --cn $CN \
    < $BCF_FILE \
    > $OUT_DIR/chrX.$CN.common_var.GWAS.txt

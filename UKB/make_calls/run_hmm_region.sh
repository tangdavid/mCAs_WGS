# set on input
# WGS_GROUP
# BATCH

set -euo pipefail

set +e 
sudo apt-get update
set -e
sudo apt --yes install parallel
sudo apt --yes install libdeflate-dev

TMP_DIR=${HOME}/tmp
mkdir -p ${TMP_DIR}

PFX=${WGS_GROUP}.batch${BATCH}
#REGION=chr14:21621904-22552132
REGION=chr1:172299011-172813287

#if [[ 1 == 0 ]]; then
echo "`date`: downloading imputed BCF..."
rm -f ${TMP_DIR}/${PFX}.imputedCommon.bcf{,.csi}
cp /mnt/project/lohdata/resources/imputedCommon/${WGS_GROUP}/${PFX}.imputedCommon.bcf{,.csi} ${TMP_DIR}
chmod u+w ${TMP_DIR}/${PFX}.imputedCommon.bcf{,.csi} 


bcftools query -l ${TMP_DIR}/${PFX}.imputedCommon.bcf \
    | parallel --joblog=/dev/stderr --colsep=' ' "NGS=WGS bash ../../tools/RAP/extract_cram_AD_RAP.sh {} ${REGION} ${TMP_DIR}/${PFX}.imputedCommon.bcf | cut -f-6 > ${TMP_DIR}/{}.${REGION}.ADs.txt"

bcftools view -r ${REGION} ${TMP_DIR}/${PFX}.imputedCommon.bcf \
    | ../../bin/mCAs_WGS prepare-bcf \
        --clear-ADs \
        --write-ADs \
        --het-format-str "${TMP_DIR}/%s.${REGION}.ADs.txt" \
        > ${TMP_DIR}/${PFX}.${REGION}.bcf 
#fi
bcftools index -f ${TMP_DIR}/${PFX}.${REGION}.bcf
rm -f ${TMP_DIR}/*.${REGION}.ADs.txt

echo "`date`: downloading CNV masks..."
zcat /mnt/project/lohdata/resources/WGS_depth/${WGS_GROUP}/${PFX}.CNVs.txt.gz \
    | grep -v MEDIAN \
    | cut -f1-5 \
    | cat - <(awk '{print "BLACKLIST\t"$0"\tMASK"}' /mnt/project/lohdata/resources/WGS_depth/lc_sv_cnv_masks.hg38.bed) \
    > ${TMP_DIR}/${PFX}.CNVs.txt

echo "`date`: downloading ref-bias..."
cp /mnt/project/lohdata/david/mCAs_WGS/cramADs_WGS/10k.WGS.ref_bias.txt ${TMP_DIR}
chmod u+w ${TMP_DIR}/10k.WGS.ref_bias.txt

../../bin/mCAs_WGS make-calls \
    --ref-bias ${TMP_DIR}/10k.WGS.ref_bias.txt \
    --cnv-mask ${TMP_DIR}/${PFX}.CNVs.txt \
    --calls-out ${OUT_DIR}/${PFX}.${REGION}.calls.txt \
    --bcf ${TMP_DIR}/${PFX}.${REGION}.bcf \
    --bcf-out /dev/null \
    --region ${REGION} \
    --phred-minor 25 \
    --phred-major 25 \
    --phred-telomere 25

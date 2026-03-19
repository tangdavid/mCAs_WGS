# set on input
# WGS_GROUP
# BATCH
# OUT_DIR
# CHR
# NGS

set -euo pipefail

set +e 
sudo apt-get update
set -e
sudo apt --yes install parallel
sudo apt --yes install libdeflate-dev

PFX=${WGS_GROUP}.batch${BATCH}
TMP_DIR=$HOME/tmp

mkdir -p $TMP_DIR

rm -f ${TMP_DIR}/${PFX}.imputedCommon.bcf{,.csi} 
cp  /mnt/project/lohdata/resources/imputedCommon/${WGS_GROUP}/${PFX}.imputedCommon.bcf{,.csi} ${TMP_DIR}

rm -f ${TMP_DIR}/GRCh38_full_analysis_set_plus_decoy_hla.fa{,.fai}
cp /mnt/project/lohdata/resources/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa{,.fai} ${TMP_DIR}

set +e
SAMPLES_FILE=${TMP_DIR}/${PFX}.samples
bcftools query -l ${TMP_DIR}/${PFX}.imputedCommon.bcf > ${SAMPLES_FILE}
cat ${SAMPLES_FILE} \
    | parallel --joblog /dev/stderr --retries 4 \
        "TMP_DIR=$TMP_DIR NGS=$NGS bash ../../tools/RAP/extract_cram_AD_RAP.sh {1} ${CHR} ${TMP_DIR}/${PFX}.imputedCommon.bcf > ${TMP_DIR}/${CHR}.{1}.ADs.txt"
set -e 

for SAMPLE in `cat ${SAMPLES_FILE}`; do
    if [[ ! -f  ${TMP_DIR}/${CHR}.${SAMPLE}.ADs.txt ]]; then 
        echo "Missing ADs for sample: " ${SAMPLE} > /dev/stderr
        continue
    fi
    awk -v OFS='\t' '{print $1,$2,$3,$4,$5,$6}' ${TMP_DIR}/${CHR}.${SAMPLE}.ADs.txt \
        > ${TMP_DIR}/${CHR}.${SAMPLE}.ADs.corrected.txt

    awk -v OFS='\t' '{print $1,$2,$3,$4,$7,$8}' ${TMP_DIR}/${CHR}.${SAMPLE}.ADs.txt \
        > ${TMP_DIR}/${CHR}.${SAMPLE}.ADs.raw.txt
done

bcftools view -r ${CHR} -Ob ${TMP_DIR}/${PFX}.imputedCommon.bcf \
    | ../../bin/mCAs_WGS prepare-bcf \
        --clear-ADs \
        --write-ADs \
        --het-format-str ${TMP_DIR}/${CHR}.%s.ADs.corrected.txt \
        > ${OUT_DIR}/${PFX}.${CHR}.correctedAD.bcf 

bcftools view -r ${CHR} -Ob ${TMP_DIR}/${PFX}.imputedCommon.bcf \
    | ../../bin/mCAs_WGS prepare-bcf \
        --clear-ADs \
        --write-ADs \
        --het-format-str ${TMP_DIR}/${CHR}.%s.ADs.raw.txt \
        > ${OUT_DIR}/${PFX}.${CHR}.rawAD.bcf 

# set on input
# WGS_GROUP
# BATCH
# OUT_DIR
# CHR
# SEX (0 for F and 1 for M)
# NGS

set -euo pipefail

set +e 
sudo apt-get update
set -e
sudo apt --yes install parallel
sudo apt --yes install libdeflate-dev

PFX=${WGS_GROUP}.batch${BATCH}
TMP_DIR=$HOME/tmp

if [[ $SEX == 1 ]]; then 
    echo "extracting chrXY ADs for M"
    CHR_NAME=chrXY
    REGION=chrX:10001-2781479,chrX:155701383-156030895
elif [[ $SEX == 0 ]]; then 
    echo "extracting chrX ADs for F"
    CHR_NAME=chrX
    REGION=chrX
else
    echo "Error: SEX has to be set as either 1 or 0"
    exit 1
fi

mkdir -p $TMP_DIR

#rm -f ${TMP_DIR}/${PFX}.imputedCommon.bcf{,.csi} 
#cp  /mnt/project/lohdata/resources/imputedCommon/${WGS_GROUP}/${PFX}.imputedCommon.bcf{,.csi} ${TMP_DIR}
bcftools view -r $REGION /mnt/project/lohdata/resources/imputedCommon/${WGS_GROUP}/${PFX}.imputedCommon.bcf \
    -Ob -o ${TMP_DIR}/${PFX}.$CHR_NAME.imputedCommon.bcf

bcftools index -f ${TMP_DIR}/${PFX}.$CHR_NAME.imputedCommon.bcf

rm -f ${TMP_DIR}/GRCh38_full_analysis_set_plus_decoy_hla.fa{,.fai}
cp /mnt/project/lohdata/resources/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa{,.fai} ${TMP_DIR}


set +e
SAMPLES_FILE=${TMP_DIR}/${PFX}.samples
grep -f \
    <(cat /mnt/project/lohdata/david/mCAs_WGS/sample_data/ID.genetic_sex.txt | awk -v SEX=$SEX '$2==SEX {print $1}' ) \
    <(bcftools query -l ${TMP_DIR}/${PFX}.$CHR_NAME.imputedCommon.bcf) \
    > ${SAMPLES_FILE}



cat ${SAMPLES_FILE} \
    | parallel --joblog /dev/stderr --retries 4 \
        "TMP_DIR=$TMP_DIR NGS=$NGS bash ../../tools/RAP/extract_cram_AD_RAP.sh {1} chrX ${TMP_DIR}/${PFX}.$CHR_NAME.imputedCommon.bcf > ${TMP_DIR}/${CHR_NAME}.{1}.ADs.txt"
set -e 

for SAMPLE in `cat ${SAMPLES_FILE}`; do
    if [[ ! -f  ${TMP_DIR}/${CHR_NAME}.${SAMPLE}.ADs.txt ]]; then 
        echo "Missing ADs for sample: " ${SAMPLE} > /dev/stderr
        continue
    fi
    awk -v OFS='\t' '{print $1,$2,$3,$4,$5,$6}' ${TMP_DIR}/${CHR_NAME}.${SAMPLE}.ADs.txt \
        > ${TMP_DIR}/${CHR_NAME}.${SAMPLE}.ADs.corrected.txt
done

bcftools view -r chrX -S $SAMPLES_FILE -Ob ${TMP_DIR}/${PFX}.$CHR_NAME.imputedCommon.bcf \
    | ../../bin/mCAs_WGS prepare-bcf \
        --clear-ADs \
        --write-ADs \
        --het-format-str ${TMP_DIR}/${CHR_NAME}.%s.ADs.corrected.txt \
        > ${OUT_DIR}/${PFX}.${CHR_NAME}.correctedAD.bcf 
bcftools index -f  ${OUT_DIR}/${PFX}.${CHR_NAME}.correctedAD.bcf

#!/bin/bash

# set on input: 
# CHR: chromosome number
# TMP_DIR: output directory and directory containing IDs
# WGS_GROUP: WGS1, WGS2, or WGS3 (for oos imputation purposes)

# output:
# $TMP_DIR/chr${CHR}.imputed.bcf*

set -euo pipefail


SNP_ARRAY_SCAFFOLD="/mnt/project/lohdata/david/mCAs_WGS/scaffold/chr${CHR}.GRCh38.40709.AC.bcf"


bcftools view $SNP_ARRAY_SCAFFOLD \
    -S $TMP_DIR/IDs.txt \
    -Ob -o $TMP_DIR/chr${CHR}.target.bcf \
    && bcftools index -f $TMP_DIR/chr${CHR}.target.bcf


REF_PANEL_NAME=`case $WGS_GROUP in
    "WGS3")
        echo "WGS1_WGS2";;
    *) echo "WGS3";;
esac`

REF="/mnt/project/lohdata/david/mCAs_WGS/imputation_reference/${REF_PANEL_NAME}/chr${CHR}.phased.${REF_PANEL_NAME}.xcf.bcf"

TARGET="$TMP_DIR/chr${CHR}.target.bcf"

MAP="/mnt/project/lohdata/david/mCAs_WGS/GRCh38_supp_files/genetic_map.chr${CHR}.txt"
CHUNKS="/mnt/project/lohdata/david/mCAs_WGS/imputation_reference/chunks_chr${CHR}.txt"

CONCAT_FILE="$TMP_DIR/chr${CHR}.concat.txt"

touch $CONCAT_FILE && rm $CONCAT_FILE
OUT_PFX=$TMP_DIR/chr${CHR}.imputed

while read CHUNK_LINE; do
    read -r CHUNK CHR_NAME BUFFER REGION <<< "$CHUNK_LINE"
    ../../bin/impute5_static \
        --h ${REF} \
        --g ${TARGET} \
        --r ${REGION} \
        --buffer-region ${BUFFER} \
        --o ${OUT_PFX}.${CHUNK}.bcf \
        --l ${OUT_PFX}.${CHUNK}.log \
        --m ${MAP} \
        --no-out-gp-field \
        --no-out-ds-field \
        --no-out-index \
        --threads `nproc`

    echo ${OUT_PFX}.${CHUNK}.bcf >> $CONCAT_FILE
done < $CHUNKS

bcftools concat -n -f $CONCAT_FILE \
    -Ob -o ${OUT_PFX}.bcf

rm ${OUT_PFX}.*.bcf
rm $TMP_DIR/chr${CHR}.target.bcf*


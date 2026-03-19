# set on input
# WGS_GROUP
# BATCH
# NGS
# OUT_DIR

set -euo pipefail
set +e 
sudo apt-get update
set -e
sudo apt --yes install libdeflate-dev

TMP_DIR=${HOME}/tmp
mkdir -p ${TMP_DIR}

PFX=${WGS_GROUP}.batch${BATCH}

WGS_GROUP=${WGS_GROUP} BATCH=${BATCH} NGS=${NGS} TMP_DIR=${TMP_DIR} bash collate_extractedADs.sh

zcat /mnt/project/lohdata/resources/WGS_depth/${WGS_GROUP}/${PFX}.CNVs.txt.gz \
    | grep -v MEDIAN \
    | cut -f1-5 \
    > ${TMP_DIR}/${PFX}.CNVs.txt

../../bin/mCAs_WGS prepare-bcf \
    --cnv-mask ${TMP_DIR}/${PFX}.CNVs.txt \
    --ref-bias-out ${OUT_DIR}/${PFX}.${NGS}.ref_bias.txt \
    < ${TMP_DIR}/${PFX}.rawAD.bcf \
    > /dev/null

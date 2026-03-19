# set on input
# WGS_GROUP
# BATCH
# OUT_DIR

set -euo pipefail

set +e 
sudo apt-get update
set -e
sudo apt --yes install parallel
PFX=${WGS_GROUP}.batch${BATCH}

TMP_DIR=${HOME}/tmp
mkdir -p ${TMP_DIR}

CALL_DIR=/mnt/project/lohdata/david/mCAs_WGS/calls

tar -xvf ${CALL_DIR}/${WGS_GROUP}/${PFX}.allelic_shift.tar -C ${TMP_DIR}
cat ${CALL_DIR}/${WGS_GROUP}/${PFX}.calls.txt > ${TMP_DIR}/${PFX}.calls.txt


for CHR in chr{1..22}; do
    echo "`date`: validating calls on ${CHR}..."
    egrep -w "ID|${CHR}" ${CALL_DIR}/${WGS_GROUP}/${PFX}.calls.txt > ${TMP_DIR}/${PFX}.${CHR}.calls.txt 
    TMP_DIR=${TMP_DIR} bash validate_calls_WES.sh \
        ${TMP_DIR}/${PFX}.${CHR}.calls.txt \
        ${TMP_DIR}/${PFX}.${CHR}.allelic_shift.bcf \
        ${TMP_DIR}/${PFX}.${CHR}.validation.WES.txt 
done

echo -e "ID\tchr\tbpStart\tbpEnd\toverRepWES\tunderRepWES\tzscore\tbdev\tbdevSE\texpectedValRate" \
    > ${OUT_DIR}/${PFX}.validation.WES.txt
cat ${TMP_DIR}/${PFX}.chr{1..22}.validation.WES.txt >> ${OUT_DIR}/${PFX}.validation.WES.txt

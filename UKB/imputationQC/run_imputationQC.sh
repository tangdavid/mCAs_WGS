# set on input: $WGS_GROUP $BATCH (smallest: WGS2 58 = 115 samples)

set -euo pipefail

set +e # OK if apt-get update fails
sudo apt-get update
set -e
sudo apt --yes install parallel
sudo apt --yes install libdeflate-dev # might be unnecessary since libdeflate0 already seems to be installed; on the other hand, doesn't hurt to make sure version is up to date

unset DX_WORKSPACE_ID # to allow dx make_download_url

PFX=${WGS_GROUP}.batch${BATCH}

TMP_DIR=$HOME/tmp
mkdir -p $TMP_DIR

ID_FILE="/mnt/project/lohdata/resources/WGS_depth/sample_batches/ID_40709.WGS_group.batch.txt"

# create sample list
awk -v WGS_GROUP=$WGS_GROUP -v BATCH=$BATCH '$2==WGS_GROUP && $3==BATCH {print $1}' \
    $ID_FILE \
    > $TMP_DIR/IDs.txt

CONCAT_FILE="$TMP_DIR/concat.txt"
touch $CONCAT_FILE && rm $CONCAT_FILE

for CHR in {1..22} X_PAR1 X; do
    CHR=$CHR TMP_DIR=$TMP_DIR bash impute_phase_chr.sh 
    echo "${TMP_DIR}/chr${CHR}.imputed.bcf" >> $CONCAT_FILE
done

set +e # OK if some jobs fail
cat $TMP_DIR/IDs.txt \
    | parallel --joblog /dev/stderr --retries 4 "TMP_DIR=$TMP_DIR bash ../../tools/RAP/extract_DRAGEN_VCF.sh {} > ${TMP_DIR}/{}.hets.txt"
set -e

bcftools concat -n -f $CONCAT_FILE -Ob -o ${TMP_DIR}/${PFX}.concat.bcf

echo "performing QC on imputed output using WGS genotypes..."
../../bin/mCAs_WGS prepare-bcf \
    --clear-ADs \
    --geno-QC \
    --het-format-str "${TMP_DIR}/%s.hets.txt" \
    < ${TMP_DIR}/${PFX}.concat.bcf \
    > ${OUT_DIR}/${PFX}.imputedCommon.bcf

bcftools index -f ${OUT_DIR}/${WGS_GROUP}.batch${BATCH}.imputedCommon.bcf


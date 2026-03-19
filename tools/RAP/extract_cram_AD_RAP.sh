# set on input: 
# NGS
# TMP_DIR

set -euo pipefail

ID=$1
REGION=$2
BCF_FILE=$3
CHR=`echo $REGION | cut -f1 -d':'`

REF=${TMP_DIR}/GRCh38_full_analysis_set_plus_decoy_hla.fa

DX_PROJECT="WGS_500K"
if [[ $NGS == "WGS" ]]; then
    CRAM_DIR="${DX_PROJECT}:/Bulk/DRAGEN WGS/Whole genome CRAM files (DRAGEN) [500k release]/${ID:0:2}"
    FIELD=24048
    CRAM_URL=`dx make_download_url "${CRAM_DIR}/${ID}_${FIELD}_0_0.dragen.cram"`
    CRAI_URL=`dx make_download_url "${CRAM_DIR}/${ID}_${FIELD}_0_0.dragen.cram.crai"`
elif [[ $NGS == "WES" ]]; then
    CRAM_DIR="${DX_PROJECT}:/Bulk/Exome sequences/Exome OQFE CRAM files/${ID:0:2}"
    FIELD=23143
    CRAM_URL=`dx make_download_url "${CRAM_DIR}/${ID}_${FIELD}_0_0.cram"`
    CRAI_URL=`dx make_download_url "${CRAM_DIR}/${ID}_${FIELD}_0_0.cram.crai"`
else
    echo "Error: NGS should be set to either WES or WGS" > /dev/stderr
    exit 1
fi 

bash `git rev-parse --show-toplevel`/tools/extract_cram_AD.sh $ID $REGION $BCF_FILE $REF $CRAM_URL $CRAI_URL

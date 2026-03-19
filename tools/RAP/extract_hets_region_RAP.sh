set -euo pipefail

ID=$1
REGION=$2

BATCH=${ID:0:2}
DX_PROJECT="WGS_500K"
VCF_DIR="${DX_PROJECT}:/Bulk/DRAGEN WGS/Whole genome variant call files (VCFs) (DRAGEN) [500k release]/$BATCH"

VCF_URL=`dx make_download_url "${VCF_DIR}/${ID}_24053_0_0.dragen.hard-filtered.vcf.gz"`
TBI_URL=`dx make_download_url "${VCF_DIR}/${ID}_24053_0_0.dragen.hard-filtered.vcf.gz.tbi"`

bcftools view -r $REGION "$VCF_URL"##idx##$TBI_URL -g het -m2 -M2 

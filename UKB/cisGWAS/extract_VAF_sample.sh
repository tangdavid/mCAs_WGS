# set on input:
# ID
# REGION

set -euo pipefail

DX_PROJECT="WGS_500K"
VCF_DIR="${DX_PROJECT}:/Bulk/DRAGEN WGS/Whole genome variant call files (VCFs) (DRAGEN) [500k release]/${ID:0:2}"
VCF_URL=`dx make_download_url "${VCF_DIR}/${ID}_24053_0_0.dragen.hard-filtered.vcf.gz"`
TBI_URL=`dx make_download_url "${VCF_DIR}/${ID}_24053_0_0.dragen.hard-filtered.vcf.gz.tbi"`
bcftools query -r $REGION "$VCF_URL"##idx##$TBI_URL -f '[%CHROM:%POS:%REF:%ALT\t%AD]\n' \
    | awk -v ID=$ID -v OFS='\t' '{print ID,$0}' \
    | tr ',' '\t' 

rm -f ${ID}_24053_0_0.dragen.hard-filtered.vcf.gz.tbi

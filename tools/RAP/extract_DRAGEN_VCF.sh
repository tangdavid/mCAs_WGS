set -euo pipefail

ID=$1

BATCH=${ID:0:2}
VCF_NAME="${ID}_24053_0_0.dragen.hard-filtered.vcf.gz"
DX_PROJECT="WGS_500K"
VCF_DIR="${DX_PROJECT}:/Bulk/DRAGEN WGS/Whole genome variant call files (VCFs) (DRAGEN) [500k release]"

`git rev-parse --show-toplevel`/bin/dragen_vcf_extract_hetADs `dx make_download_url "${VCF_DIR}/${BATCH}/${VCF_NAME}"` \
    | awk '$1 ~/chr([1-9]|1[0-9]|2[0-2]|X)$/' 


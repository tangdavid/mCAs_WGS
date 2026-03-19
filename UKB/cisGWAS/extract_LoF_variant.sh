# set on input 
# REGION
# GENE
# MASK
# BURDEN_MASKS_DIR
# ID

set -euo pipefail

CHR=`echo $REGION | cut -f1 -d:`

DX_PROJECT="WGS_500K"
VCF_DIR="${DX_PROJECT}:/Bulk/DRAGEN WGS/Whole genome variant call files (VCFs) (DRAGEN) [500k release]/${ID:0:2}"
VCF_URL=`dx make_download_url "${VCF_DIR}/${ID}_24053_0_0.dragen.hard-filtered.vcf.gz"`
TBI_URL=`dx make_download_url "${VCF_DIR}/${ID}_24053_0_0.dragen.hard-filtered.vcf.gz.tbi"`

PREFIX="$ID."`echo $REGION | tr ':-' '_'`
bcftools view -r $REGION "$VCF_URL"##idx##$TBI_URL -g het -m2 -M2 \
    | bcftools query -f '[%CHROM\t%POS\t%REF\t%ALT\n]' \
    | awk -v MASK=$MASK -v GENE=$GENE -v ID=$ID -v OFS='\t' '
ARGIND==1 && $1==MASK {
    split($2, a, ",");
    for(i in a) {
        mask_elements[a[i]] = 1;
    }
}
ARGIND==2 && $2==GENE && ($3 in mask_elements)  {
    save["chr"$1]=1
}
ARGIND==3 && ($1":"$2":"$3":"$4 in save) {
    print ID,$0
}' \
    $BURDEN_MASKS_DIR/$CHR.gnomAD.LoF.missense.CNV.masks \
    <(zcat $BURDEN_MASKS_DIR/$CHR.gnomAD.LoF.missense.CNV.annotations.txt.gz | grep -w $GENE ) \
    - 

rm -f ${ID}_24053_0_0.dragen.hard-filtered.vcf.gz.tbi

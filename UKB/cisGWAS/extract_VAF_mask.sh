# require gnu-parallel to be installed
# set on input
# CHR
# GENE
# MASK
# AF

set -euo pipefail

TMP_DIR=$HOME/tmp
mkdir -p $TMP_DIR

GENE_MASK_AF=$GENE.$MASK.$AF
RARE_VAR_VCF=${RARE_VAR_VCF:-/mnt/project/lohdata/resources/burden_masks/vcfs/$CHR.LoF.scaffold.phased.AS.bcf}
BURDEN_MASKS_DIR=${BURDEN_MASKS_DIR:-/mnt/project/lohdata/resources/burden_masks/mask_definitions}
PREFIX=$GENE_MASK_AF

# extract a region that contains all the protein coding variants for the gene and nearby common variants
REGION=`zcat $BURDEN_MASKS_DIR/gnomAD.protein_coding.high_moderate_vep.$CHR.txt.gz \
    | awk -v GENE=$GENE '
BEGIN {
    minpos=1e9;
    maxpos=0;
}
$10==GENE {
    chrom=$1;
    minpos=($2<minpos) ? $2: minpos; 
    maxpos=($2>maxpos) ? $2 : maxpos
} 
END {
    print chrom":"minpos"-"maxpos
}'`

# extract LoF variant for all individuals
echo 'Extracting LoF variant...'
bcftools query -r $REGION -f '[%SAMPLE\t%ID\t%GT\n]' $RARE_VAR_VCF \
    | grep -w $GENE_MASK_AF \
    | awk 'ARGIND==1 {keep[$1]=1} ARGIND==2 && ($1 in keep) && ($3=="1|0" || $3=="0|1" || $3=="0/1" || $3=="1/0") {print $1}' \
        /mnt/project/lohdata/david/mCAs_WGS/sample_data/WGS_500k.included_IDs.txt \
        - \
    | parallel --joblog=/dev/stderr --halt-on-error 2 "REGION=$REGION GENE=$GENE MASK=$MASK BURDEN_MASKS_DIR=$BURDEN_MASKS_DIR ID={} bash extract_LoF_variant.sh" \
    > $TMP_DIR/$PREFIX.LoF_variant.txt 

echo 'Extracting VAF...'
cat $TMP_DIR/$PREFIX.LoF_variant.txt \
    | awk -v OFS='\t' '{print $1,$2":"$3"-"$3}' \
    | parallel --joblog=/dev/stderr --halt-on-error 2 --colsep='\t' "ID={1} REGION={2} bash extract_VAF_sample.sh" \
    > $TMP_DIR/$PREFIX.VAF.txt

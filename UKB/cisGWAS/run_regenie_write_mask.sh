# set on input:
# CHR_NUM
# OUT_DIR

set -euo pipefail

TMP_DIR=$HOME/tmp
mkdir -p $TMP_DIR

#BGEN_PREFIX="/mnt/project/Bulk/Exome sequences/Population level exome OQFE variants, BGEN format - final release/ukb23159_c${CHR_NUM}_b0_v1"
BGEN_PREFIX="/mnt/project/Bulk/DRAGEN WGS/DRAGEN population level WGS variants, BGEN format [500k release]/ukb24309_c${CHR_NUM}_b0_v1"
MASKS_DIR="/mnt/project/lohdata/resources/burden_masks/mask_definitions"
WES_CNV_BGEN_PREFIX="/mnt/project/lohdata/resources/WES_CNVs/chr$CHR_NUM.gene.v3"
SCAFFOLD_VCF="/mnt/project/lohdata/david/mCAs_WGS/scaffold/chr$CHR_NUM.GRCh38.40709.AC.bcf"
WITHDRAWN="/mnt/project/lohdata/resources/BOLT-LMM/w40709_CURRENT.FID_IID.txt"

# WGS bgens are stored ref-last
./extractBgenVariants \
    <(zcat $MASKS_DIR/gnomAD.protein_coding.high_moderate_vep.chr$CHR_NUM.txt.gz \
        | awk -v OFS='\t' '
        {   
            chr_num=substr($1,4);
            if(chr_num=="X" && $2<=2781479) chr_num="PAR1";
            else if(chr_num=="X" && $2>=155701383) chr_num="PAR2";
            print chr_num,$2,$5,$4
        }
        '
    ) \
    "$BGEN_PREFIX.bgen" \
    $TMP_DIR/chr$CHR_NUM.protein_coding.bgen
bgenix -g $TMP_DIR/chr$CHR_NUM.protein_coding.bgen -index -clobber

# Output file prefixes
PROTEIN_CODING_PREFIX=$TMP_DIR/chr$CHR_NUM.protein_coding
CNV_PREFIX=$TMP_DIR/chr$CHR_NUM.CNV
MERGED_PREFIX=$TMP_DIR/chr$CHR_NUM.protein_coding.CNV
BURDEN_MASKS_DIR=$OUT_DIR/mask_definitions

plink2 --bgen $TMP_DIR/chr$CHR_NUM.protein_coding.bgen ref-last \
    --sample "$BGEN_PREFIX.sample" \
    --rm-dup force-first \
    --hard-call-threshold 0.25 \
    --remove $WITHDRAWN \
    --set-all-var-ids @:#:\$r:\$a \
    --new-id-max-allele-len 1000 missing \
    --geno 0.01 \
    --mac 1 \
    --make-bed --out $PROTEIN_CODING_PREFIX

if [[ $CHR_NUM == "X" ]]; then
    sed -i 's/PAR1\|PAR2/X/g' $PROTEIN_CODING_PREFIX.bim 
    MERGED_PREFIX=$PROTEIN_CODING_PREFIX
else 
    plink2 --bgen $WES_CNV_BGEN_PREFIX.bgen ref-first \
        --sample $WES_CNV_BGEN_PREFIX.sample \
        --hard-call-threshold 0.25 \
        --remove $WITHDRAWN \
        --make-bed --out $CNV_PREFIX

    plink --bfile $CNV_PREFIX \
        --bmerge $PROTEIN_CODING_PREFIX \
        --make-bed --out $MERGED_PREFIX
fi

mkdir -p $BURDEN_MASKS_DIR
CHR=chr$CHR_NUM BFILE=$PROTEIN_CODING_PREFIX OUT_DIR=$BURDEN_MASKS_DIR bash make_burden_masks.sh

PHENO_FILE=$TMP_DIR/pheno_file.txt
awk 'NR!=2' "$BGEN_PREFIX.sample" | sed 's/ID_1 ID_2/FID IID/' > $PHENO_FILE
BURDEN_MASKS_DEF=$BURDEN_MASKS_DIR/chr$CHR_NUM.gnomAD.LoF.missense.CNV

BURDEN_MASKS_PREFIX=$TMP_DIR/chr$CHR_NUM.LoF.missense.CNV

regenie \
    --step 2 \
    --ignore-pred \
    --bed $MERGED_PREFIX \
    --phenoFile $PHENO_FILE \
    --set-list $BURDEN_MASKS_DEF.sets.txt.gz \
    --anno-file $BURDEN_MASKS_DEF.annotations.txt.gz \
    --mask-def $BURDEN_MASKS_DEF.masks \
    --nauto 22 \
    --aaf-bins 0.01,0.001,0.0001 \
    --bsize 200 \
    --write-mask \
    --skip-test \
    --out $BURDEN_MASKS_PREFIX

# preparing bcf with burden masks
plink2 \
    --bfile ${BURDEN_MASKS_PREFIX}_masks \
    --export bcf \
    --mind 0.1 \
    --output-chr chrM \
    --no-fid \
    --out $BURDEN_MASKS_PREFIX.masks

mkdir -p $OUT_DIR/bgen
plink2 \
    --bfile ${BURDEN_MASKS_PREFIX}_masks \
    --export bgen-1.2 bits=8 \
    --out $OUT_DIR/bgen/chr$CHR_NUM.burden_masks
bgenix -g $OUT_DIR/bgen/chr$CHR_NUM.burden_masks.bgen -index -clobber

if [[ $CHR_NUM == "X" ]]; then 
    exit
fi

# preparing bcf with single variants
plink2 \
    --bfile $PROTEIN_CODING_PREFIX \
    --export bcf \
    --output-chr chrM \
    --mac 10 \
    --no-fid \
    --nonfounders \
    --extract <(zcat $MASKS_DIR/gnomAD.protein_coding.high_moderate_vep.chr$CHR_NUM.txt.gz | awk '{print substr($1,4)":"$2":"$4":"$5}') \
    --out ${PROTEIN_CODING_PREFIX}.MAC10

mkdir -p $OUT_DIR/vcfs/tmp
SCAFFOLD_SUBSET=$OUT_DIR/vcfs/tmp/chr$CHR_NUM.snp_scaffold.subset_IDs.bcf

echo "`date`: get common samples between burden masks and scaffold..."
COMMON_SAMPLES=$TMP_DIR/common_sample_IDs.txt
SAMPLE_ORDER=$TMP_DIR/sample_order.txt
awk 'ARGIND==1 {keep[$1]=1} ARGIND==2 && keep[$1]==1 {print $1}' \
    <(bcftools query -l $BURDEN_MASKS_PREFIX.masks.bcf) \
    <(bcftools query -l ${PROTEIN_CODING_PREFIX}.MAC10.bcf) \
    | awk 'ARGIND==1 {keep[$1]=1} ARGIND==2 && keep[$1]==1 {print $1}' \
    - \
    <(bcftools query -l $SCAFFOLD_VCF) \
    | sort -u \
    > $COMMON_SAMPLES

# subset scaffold VCF to common samples
echo "`date`: subset scaffold VCF to common samples..."
bcftools view \
    -i 'INFO/AC>=10000 && INFO/AC<=990000' \
    -S $COMMON_SAMPLES \
    -Ob -o $SCAFFOLD_SUBSET \
    $SCAFFOLD_VCF
bcftools index -f $SCAFFOLD_SUBSET

# reorder samples to match the order in the scaffold file
bcftools query -l $SCAFFOLD_SUBSET > $SAMPLE_ORDER

echo "`date`: reorder samples in burden masks bcf..."
bcftools view -S $SAMPLE_ORDER -Ou $BURDEN_MASKS_PREFIX.masks.bcf \
    | bcftools +fill-tags -Ob -o $BURDEN_MASKS_PREFIX.masks.reorder.bcf -- -t AC,AN
bcftools index -f $BURDEN_MASKS_PREFIX.masks.reorder.bcf

echo "`date`: reorder samples in single variants bcf..."
bcftools view -S $SAMPLE_ORDER -Ou ${PROTEIN_CODING_PREFIX}.MAC10.bcf \
    | bcftools +fill-tags -Ob -o ${PROTEIN_CODING_PREFIX}.MAC10.reorder.bcf -- -t AC,AN
bcftools index -f ${PROTEIN_CODING_PREFIX}.MAC10.reorder.bcf

# merging all bcf files for phase_rare
echo "`date`: merge burden masks and single variants with haplotype scaffold..."
bcftools concat --allow-overlaps -Ob \
    $BURDEN_MASKS_PREFIX.masks.reorder.bcf \
    ${PROTEIN_CODING_PREFIX}.MAC10.reorder.bcf \
    $SCAFFOLD_SUBSET \
    -o $OUT_DIR/vcfs/tmp/chr$CHR_NUM.phase_rare.input.bcf
bcftools index -f $OUT_DIR/vcfs/tmp/chr$CHR_NUM.phase_rare.input.bcf

echo "`date`: done!"

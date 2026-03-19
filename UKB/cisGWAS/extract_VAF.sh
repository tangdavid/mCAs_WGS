set -euo pipefail

PHASING_DIR=/mnt/project/lohdata/david/mCAs_WGS/cisGWAS/results/refined_phase

cat /mnt/project/lohdata/david/mCAs_WGS/cisGWAS/results/CN-LOH.GWAS.FDR.01.txt  \
    | cut -f3 \
    | grep -wf - <(zcat $PHASING_DIR/ID.variant.mask.read_phase.statistical_phase.chr{1..22}.txt.gz ) \
    | awk -v OFS='\t' '$2!="." {split($2,a,":"); print $1,$2,a[1]":"a[2]"-"a[2]}' \
    | tee variants.txt \
    | parallel --colsep='\t' --joblog=/dev/stderr "ID={1} REGION={3} bash extract_VAF_sample.sh" > ADs.txt 


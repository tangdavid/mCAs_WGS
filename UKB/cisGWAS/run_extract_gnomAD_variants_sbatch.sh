set -euo pipefail

for CHR in chr{1..22}
do
    MEM=1G

    OUT_DIR=/n/data1/bwh/medicine/loh/david/mCAs_WGS/cisGWAS/burden_masks
    sbatch \
        -p short \
        -t 1:0:0 \
        -c 1 \
        --mem $MEM \
        --account loh \
        -J $CHR.extract_gnomAD_variants \
        -o /n/data1/bwh/medicine/loh/david/mCAs_WGS/cisGWAS/burden_masks/run_extract_gnomAD_variants.${CHR}.sbatch.log \
        --wrap "srun /usr/bin/time -v bash run_extract_gnomAD_variants.sh $CHR $OUT_DIR"
done

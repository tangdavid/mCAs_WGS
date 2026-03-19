# set on input:
# CHR

set -euo pipefail

DIR="/lohdata/david/mCAs_WGS"

echo $CHR

dx run /lohdata/resources/hts-knife \
    -iin="$DIR"/run_pipeline.sh \
    -icmd="set -euo pipefail && CHR=$CHR /usr/bin/time -v bash run_pipeline.sh cisGWAS run_extract_gnomAD_variants.sh 2>&1 | tee $CHR.extract_gnomAD_variants.dx_run.log" \
    --instance-type mem2_ssd1_v2_x4 \
    --destination "$DIR"/cisGWAS/burden_masks/ \
    --priority low \
    --brief \
    --ignore-reuse \
    --name extract_gnomAD_variants.$CHR \
    --allow-ssh \
    --tag extract_gnomAD_variants \
    -y

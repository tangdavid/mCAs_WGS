# set on input:
# CHR_NUM

set -euo pipefail

DIR="/lohdata/david/mCAs_WGS"

echo chr$CHR_NUM

dx run swiss-army-knife \
    -iin="$DIR"/run_pipeline.sh \
    -icmd="set -euo pipefail && CHR_NUM=$CHR_NUM /usr/bin/time -v bash run_pipeline.sh cisGWAS run_regenie_write_mask.sh 2>&1 | tee chr$CHR_NUM.write_masks.dx_run.log" \
    --instance-type mem1_ssd2_v2_x4 \
    --destination /lohdata/resources/burden_masks/ \
    --priority high \
    --brief \
    --ignore-reuse \
    --name run_regenie_write_mask.chr$CHR_NUM \
    --allow-ssh \
    --tag run_regenie_write_mask \
    -y

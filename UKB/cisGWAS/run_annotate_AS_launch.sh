# set on input:
# CHR

set -euo pipefail

DIR="/lohdata/david/mCAs_WGS"


dx run /lohdata/resources/hts-knife \
    -iin="$DIR"/run_pipeline.sh \
    -icmd="set -euo pipefail && CHR=$CHR /usr/bin/time -v bash run_pipeline.sh cisGWAS annotate_AS.sh 2>&1 | tee $CHR.annotate_AS.dx_run.log" \
    --instance-type mem1_ssd1_v2_x4 \
    --destination /lohdata/resources/burden_masks/vcfs/ \
    --priority low \
    --brief \
    --ignore-reuse \
    --name annotate_AS.$CHR \
    --allow-ssh \
    --extra-args '{"executionPolicy": {"maxRestarts": 5}}' \
    --tag annotate_AS \
    -y

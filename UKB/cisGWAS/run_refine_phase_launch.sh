# set on input:
# CHR

set -euo pipefail

DIR="/lohdata/david/mCAs_WGS"

echo $CHR

dx run /lohdata/resources/hts-knife \
    -iin="$DIR"/run_pipeline.sh \
    -icmd="set -euo pipefail && CHR=$CHR /usr/bin/time -v bash run_pipeline.sh cisGWAS run_refine_phase.sh 2>&1 | tee $CHR.refine_phase.dx_run.log" \
    --instance-type mem1_ssd1_v2_x4 \
    --destination "$DIR"/cisGWAS/results/refined_phase/ \
    --priority low \
    --brief \
    --ignore-reuse \
    --name refine_phase.$CHR \
    --allow-ssh \
    --extra-args '{"executionPolicy": {"maxRestarts": 5}}' \
    --tag refine_phase \
    -y

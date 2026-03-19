set -euo pipefail

DIR="/lohdata/david/mCAs_WGS"

dx run swiss-army-knife \
    -iin="$DIR"/run_pipeline.sh \
    -icmd="set -euo pipefail && /usr/bin/time -v bash run_pipeline.sh transGWAS test_13q_CNLOH.sh 2>&1 | tee regenie.13qCNLOH.dx_run.log" \
    --instance-type mem1_ssd1_v2_x16 \
    --destination "$DIR"/del13q/13qCNLOH/ \
    --priority high \
    --brief \
    --ignore-reuse \
    --name regenie.13qCNLOH \
    --allow-ssh \
    -y

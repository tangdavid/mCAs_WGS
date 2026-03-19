# set on input:

set -euo pipefail

DIR="/lohdata/david/mCAs_WGS"

dx run swiss-army-knife \
    -iin="$DIR"/run_pipeline.sh \
    -icmd="set -euo pipefail && /usr/bin/time -v bash run_pipeline.sh transGWAS regenie_CLL_step1.sh 2>&1 | tee regenie.CLL.step1.dx_run.log" \
    --instance-type mem1_ssd1_v2_x16 \
    --destination "$DIR"/del13q/GWAS/ \
    --priority high \
    --brief \
    --ignore-reuse \
    --name regenie.CLL.step1 \
    --allow-ssh \
    -y

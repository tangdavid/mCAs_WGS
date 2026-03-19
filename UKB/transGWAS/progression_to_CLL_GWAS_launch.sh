# set on input:
# CHR

set -euo pipefail

DIR="/lohdata/david/mCAs_WGS"

dx run swiss-army-knife \
    -iin="$DIR"/run_pipeline.sh \
    -icmd="set -euo pipefail && CHR=$CHR /usr/bin/time -v bash run_pipeline.sh transGWAS progression_to_CLL_GWAS.sh 2>&1 | tee $CHR.regenie.CLL_progression.dx_run.log" \
    --instance-type mem1_ssd1_v2_x16 \
    --destination "$DIR"/del13q/CLL_progression/ \
    --priority high \
    --brief \
    --ignore-reuse \
    --name $CHR.regenie.CLL_progression \
    --allow-ssh \
    -y

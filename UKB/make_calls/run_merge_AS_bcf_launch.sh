set -euo pipefail

DIR="/lohdata/david/mCAs_WGS"
CORES=8

dx run swiss-army-knife \
    -iin="$DIR"/run_pipeline.sh \
    -icmd="set -euo pipefail && /usr/bin/time -v bash run_pipeline.sh make_calls merge_AS_bcf.sh 2>&1 | tee merge_AS_bcf.dx_run.log" \
    --instance-type mem1_ssd1_v2_x$CORES \
    --destination "$DIR"/calls/allelic_shift/ \
    --priority high \
    --brief \
    --ignore-reuse \
    --name merge_AS_bcf \
    --allow-ssh \
    -y

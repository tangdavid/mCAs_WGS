# set on input:
# WGS_GROUP
# BATCH
# NAME

set -euo pipefail

DIR="/lohdata/david/mCAs_WGS"

CORES=8

echo $WGS_GROUP $BATCH $NAME

dx run /lohdata/resources/hts-knife \
    -iin="$DIR"/run_pipeline.sh \
    -icmd="set -euo pipefail && WGS_GROUP=$WGS_GROUP BATCH=$BATCH /usr/bin/time -v bash run_pipeline.sh benchmark_mocha run_benchmark_mocha.sh 2>&1 | tee $WGS_GROUP.batch$BATCH.mocha.dx_run.log" \
    --instance-type mem1_ssd1_v2_x$CORES \
    --destination "$DIR"/benchmark_mocha/$WGS_GROUP \
    --priority high \
    --brief \
    --ignore-reuse \
    --name $NAME.$WGS_GROUP.batch$BATCH \
    --allow-ssh \
    --tag $NAME \
    --tag $NAME.$WGS_GROUP \
    -y

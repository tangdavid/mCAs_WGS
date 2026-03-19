set -euo pipefail

DIR="/lohdata/david/mCAs_WGS"

CORES=8
echo $WGS_GROUP $BATCH $NAME

dx run /lohdata/resources/hts-knife \
    -iin="$DIR"/run_pipeline.sh \
    -icmd="set -euo pipefail && WGS_GROUP=$WGS_GROUP BATCH=$BATCH /usr/bin/time -v bash run_pipeline.sh make_calls run_all_calls.sh 2>&1 | tee $WGS_GROUP.batch$BATCH.dx_run.log" \
    --instance-type mem1_ssd1_v2_x$CORES \
    --destination "$DIR"/pilot10k/ \
    --priority low \
    --brief \
    --ignore-reuse \
    --name $NAME.$WGS_GROUP.batch$BATCH \
    --allow-ssh \
    -y

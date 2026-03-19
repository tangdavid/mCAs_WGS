# set on input:
# WGS_GROUP
# BATCH
# NAME

set -euo pipefail

DIR="/lohdata/david/mCAs_WGS"

CORES=4

echo $WGS_GROUP $BATCH $NAME

dx run /lohdata/resources/hts-knife \
    -iin="$DIR"/run_pipeline.sh \
    -icmd="set -euo pipefail && WGS_GROUP=$WGS_GROUP BATCH=$BATCH /usr/bin/time -v bash run_pipeline.sh make_calls run_validate_calls.sh 2>&1 | tee $WGS_GROUP.batch$BATCH.$NAME.dx_run.log" \
    --instance-type mem1_ssd1_v2_x$CORES \
    --destination "$DIR"/calls/$WGS_GROUP \
    --priority low \
    --brief \
    --ignore-reuse \
    --name $NAME.$WGS_GROUP.batch$BATCH \
    --allow-ssh \
    --tag $NAME \
    --tag $NAME.$WGS_GROUP \
    -y

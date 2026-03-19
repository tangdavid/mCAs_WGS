set -euo pipefail

DIR="/lohdata/david/mCAs_WGS"

CORES=16

echo $WGS_GROUP $BATCH $NAME

NGS=WGS 

dx run /lohdata/resources/hts-knife \
    -iin="$DIR"/run_pipeline.sh \
    -icmd="set -euo pipefail && WGS_GROUP=$WGS_GROUP BATCH=$BATCH NGS=$NGS SEX=0 /usr/bin/time -v bash run_pipeline.sh mLOXY run_extractADs_chrX.sh 2>&1 | tee $NAME.$WGS_GROUP.batch$BATCH.dx_run.log" \
    --instance-type mem1_ssd1_v2_x$CORES \
    --destination "$DIR"/$NAME/$WGS_GROUP/ \
    --priority high \
    --brief \
    --ignore-reuse \
    --name $NAME.$WGS_GROUP.batch$BATCH \
    --tag $NAME \
    --allow-ssh \
    -y

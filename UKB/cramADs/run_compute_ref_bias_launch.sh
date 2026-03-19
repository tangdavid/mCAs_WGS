# set on input:
# WGS_GROUP
# BATCH
# NAME

set -euo pipefail

DIR="/lohdata/david/mCAs_WGS"

CORES=2
NGS=`echo $NAME | cut -f2 -d.`

if [ $NGS != "WGS" ] && [ $NGS != "WES" ]; then 
    echo "error: must specify either WGS or WES in the name"
    exit 1
fi 

echo $WGS_GROUP $BATCH $NAME

dx run /lohdata/resources/hts-knife \
    -iin="$DIR"/run_pipeline.sh \
    -icmd="set -euo pipefail && WGS_GROUP=$WGS_GROUP BATCH=$BATCH NGS=$NGS /usr/bin/time -v bash run_pipeline.sh cramADs run_compute_ref_bias.sh 2>&1 | tee $WGS_GROUP.batch$BATCH.$NAME.dx_run.log" \
    --instance-type mem1_ssd1_v2_x$CORES \
    --destination "$DIR"/cramADs_${NGS}/$WGS_GROUP \
    --priority low \
    --brief \
    --ignore-reuse \
    --name $NAME.$WGS_GROUP.batch$BATCH \
    --allow-ssh \
    --tag $NAME \
    --tag $NAME.$WGS_GROUP \
    -y

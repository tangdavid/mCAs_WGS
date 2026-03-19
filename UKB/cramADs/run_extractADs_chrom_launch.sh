set -euo pipefail

DIR="/lohdata/david/mCAs_WGS"

CORES=8
NGS=`echo $NAME | cut -f2 -d.`

if [ $NGS != "WGS" ] && [ $NGS != "WES" ]; then 
    echo "error: must specify either WGS or WES in the name"
    exit 1
fi 

echo $WGS_GROUP $BATCH $NAME

for CHR in chr{{1..22},X}
do
    if [ `egrep 'done|running|runnable' ../$NAME.job_status.txt | grep -c "$NAME.$WGS_GROUP.batch$BATCH.$CHR"` != 0 ]; then continue; fi
    dx run /lohdata/resources/hts-knife \
        -iin="$DIR"/run_pipeline.sh \
        -icmd="set -euo pipefail && WGS_GROUP=$WGS_GROUP BATCH=$BATCH NGS=$NGS CHR=$CHR /usr/bin/time -v bash run_pipeline.sh cramADs run_extractADs_chrom.sh 2>&1 | tee $WGS_GROUP.batch$BATCH.cramADs.$CHR.dx_run.log" \
        --instance-type mem1_ssd1_v2_x$CORES \
        --destination "$DIR"/cramADs_${NGS}/$WGS_GROUP/ \
        --priority low \
        --brief \
        --ignore-reuse \
        --name $NAME.$WGS_GROUP.batch$BATCH.$CHR \
        --tag $NAME \
        --allow-ssh \
        -y
done


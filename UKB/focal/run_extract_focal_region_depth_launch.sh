set -euo pipefail

DIR="/lohdata/david/mCAs_WGS"

grep -w $FOCAL_REGION focal_region_list.txt || exit
for FIRST_BATCH in {0..49}; do
    echo "FIRST_BATCH: $FIRST_BATCH"
    dx run /lohdata/resources/hts-knife \
        -iin="$DIR"/run_pipeline.sh \
        -icmd="set -euo pipefail && FIRST_BATCH=$FIRST_BATCH FOCAL_REGION=$FOCAL_REGION /usr/bin/time -v bash run_pipeline.sh focal run_extract_focal_region_depth.sh 2>&1 | tee $FIRST_BATCH.dx_run.log" \
        --instance-type mem2_ssd1_v2_x2 \
        --destination "$DIR"/$FOCAL_REGION/depth/batches/ \
        --priority low \
        --brief \
        --ignore-reuse \
        --name $FOCAL_REGION.depth.$FIRST_BATCH \
        --allow-ssh \
        --tag $FOCAL_REGION.depth \
        -y
done

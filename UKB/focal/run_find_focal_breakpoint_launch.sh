set -euo pipefail

DIR="/lohdata/david/mCAs_WGS"
FOCAL_REGION=${1:-del13q}
MAX_LENGTH=${MAX_LENGTH:-5000}
MIN_LENGTH=${MIN_LENGTH:-100}

grep -w $FOCAL_REGION focal_region_list.txt || exit

for FIRST_BATCH in {0..49}; do
    echo "FIRST_BATCH: $FIRST_BATCH"
    dx run /lohdata/resources/hts-knife \
        -iin="$DIR"/run_pipeline.sh \
        -icmd="set -euo pipefail && FIRST_BATCH=$FIRST_BATCH FOCAL_REGION=$FOCAL_REGION MAX_LENGTH=$MAX_LENGTH MIN_LENGTH=$MIN_LENGTH /usr/bin/time -v bash run_pipeline.sh focal run_find_focal_breakpoint.sh 2>&1 | tee $FOCAL_REGION.$FIRST_BATCH.dx_run.log" \
        --instance-type mem2_ssd1_v2_x2 \
        --destination "$DIR"/$FOCAL_REGION/breakpoints/batches/ \
        --priority high \
        --brief \
        --ignore-reuse \
        --name $FOCAL_REGION.$FIRST_BATCH.breakpoints \
        --allow-ssh \
        --tag $FOCAL_REGION.breakpoints \
        -y
done

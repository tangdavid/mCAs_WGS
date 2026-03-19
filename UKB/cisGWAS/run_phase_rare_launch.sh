# set on input:
# CHR

set -euo pipefail

DIR="/lohdata/david/mCAs_WGS"

dx ls "/lohdata/resources/burden_masks/vcfs/tmp/" | grep -q $CHR.snp_scaffold.subset_IDs.bcf || exit

echo $CHR
for BATCH in `dx ls "$DIR/imputation_reference/small_chunks_6cM/jobs" | grep -w chunks_${CHR} | cut -f2 -d.`; do 
    echo $CHR $BATCH
    dx run  /lohdata/resources/hts-knife \
        -iin="$DIR"/run_pipeline.sh \
        -icmd="set -euo pipefail && CHR=$CHR BATCH=$BATCH /usr/bin/time -v bash run_pipeline.sh cisGWAS run_phase_rare.sh 2>&1 | tee $CHR.b$BATCH.phase_rare.dx_run.log" \
        --instance-type mem1_ssd1_v2_x16 \
        --destination /lohdata/resources/burden_masks/vcfs/ \
        --priority low \
        --brief \
        --ignore-reuse \
        --name phase_rare.$CHR.b$BATCH \
        --extra-args '{"executionPolicy": {"maxRestarts": 5}}' \
        --allow-ssh \
        --tag phase_rare \
        -y
done

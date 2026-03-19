set -euo pipefail

set +e 
sudo apt-get update
set -e
sudo apt --yes install parallel

export TMP_DIR=$HOME/tmp
mkdir -p $TMP_DIR
rm -f $TMP_DIR/WGS{1..6}.mLOX.files.txt

source /mnt/project/lohdata/david/install_mocha.sh

function extract_calls() {
    # set on input: WGS_GROUP, BATCH
    PFX=$WGS_GROUP.batch$BATCH
    MLOX_DIR=/mnt/project/lohdata/david/mCAs_WGS/mLOX_new_phase/$WGS_GROUP
    bcftools view $MLOX_DIR/$PFX.mLOX.allelic_shift.bcf \
        -S <(zcat $MLOX_DIR/$PFX.mLOX.calls.txt.gz | awk 'NR>1' | cut -f1 | sort -u) \
        -Ou \
        | bcftools +extendFMT --format AS --phase --dist 500000 -Ou \
        > $TMP_DIR/$PFX.mLOX.calls.bcf
}

function merge_WGS_GROUP() {
    # set on input: WGS_GROUP
    bcftools merge --force-single -m none -l $TMP_DIR/$WGS_GROUP.mLOX.files.txt -Ob --no-index > $TMP_DIR/$WGS_GROUP.mLOX.merged.bcf
    cat $TMP_DIR/$WGS_GROUP.mLOX.files.txt | xargs -I{} rm {}
    rm -f $TMP_DIR/$WGS_GROUP.mLOX.files.txt
}

export -f extract_calls 
export -f merge_WGS_GROUP

BATCH_INFO=/mnt/project/lohdata/resources/WGS_depth/sample_batches/WGS_group.batch.Nsamples.txt
for WGS_GROUP in WGS{1..6}; do
    grep $WGS_GROUP $BATCH_INFO | awk -v TMP_DIR=$TMP_DIR '!($1=="WGS5" && $2==127) {print TMP_DIR"/"$1".batch"$2".mLOX.calls.bcf"}' > $TMP_DIR/$WGS_GROUP.mLOX.files.txt
    grep $WGS_GROUP $BATCH_INFO \
        | awk '!($1=="WGS5" && $2==127)' \
        | parallel --colsep='\t' --joblog=/dev/stderr WGS_GROUP={1} BATCH={2} extract_calls 

    WGS_GROUP=$WGS_GROUP merge_WGS_GROUP
done

bcftools merge -m none $TMP_DIR/WGS{1..6}.mLOX.merged.bcf -Ob --no-index > $OUT_DIR/WGS_500k.mLOX.allelic_shift.bcf
bcftools index -f $OUT_DIR/WGS_500k.mLOX.allelic_shift.bcf
rm $TMP_DIR/WGS{1..6}.mLOX.merged.bcf


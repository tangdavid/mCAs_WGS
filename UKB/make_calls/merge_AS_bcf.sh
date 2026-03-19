set -euo pipefail

set +e 
sudo apt-get update
set -e
sudo apt --yes install parallel

export TMP_DIR=$HOME/tmp
mkdir -p $TMP_DIR
rm -f $TMP_DIR/WGS{1..6}.chr{1..22}.files.txt

source /mnt/project/lohdata/david/install_mocha.sh

function extract_calls() {
    CHR=$1
    if ! zgrep -wq $CHR /mnt/project/lohdata/david/mCAs_WGS/calls/$WGS_GROUP/$PFX.calls.txt.gz; then
        bcftools view -Ou -G $TMP_DIR/$PFX.$CHR.allelic_shift.bcf >  $TMP_DIR/$PFX.$CHR.calls.bcf
    else 
        zcat /mnt/project/lohdata/david/mCAs_WGS/calls/$WGS_GROUP/$PFX.calls.txt.gz | grep -w $CHR | cut -f1 | sort -u > $TMP_DIR/$PFX.$CHR.samples
        zgrep -wq $CHR /mnt/project/lohdata/david/mCAs_WGS/calls/$WGS_GROUP/$PFX.calls.txt.gz || rm $TMP_DIR/$PFX.$CHR.allelic_shift.bcf{,.csi} 
        bcftools view $TMP_DIR/$PFX.$CHR.allelic_shift.bcf \
            -S $TMP_DIR/$PFX.$CHR.samples -Ou \
            | bcftools +extendFMT --format AS --phase --dist 500000 -Ou \
            > $TMP_DIR/$PFX.$CHR.calls.bcf
        rm $TMP_DIR/$PFX.$CHR.samples
    fi
    echo $TMP_DIR/$PFX.$CHR.calls.bcf >> $TMP_DIR/$WGS_GROUP.$CHR.files.txt 
    rm $TMP_DIR/$PFX.$CHR.allelic_shift.bcf{,.csi}
}

function merge_WGS_GROUP() {
    echo chr{1..22} \
        | tr ' ' '\n' \
        | parallel --halt-on-error 2 --joblog=/dev/stderr "bcftools merge --force-single -m none -l $TMP_DIR/$WGS_GROUP.{}.files.txt -Ob --no-index > $TMP_DIR/$WGS_GROUP.{}.merged.bcf"
    cat $TMP_DIR/$WGS_GROUP.chr{1..22}.files.txt | xargs -I{} rm {}
    rm -f $TMP_DIR/$WGS_GROUP.chr{1..22}.files.txt
}

export -f extract_calls
export -f merge_WGS_GROUP

PREV_GROUP=WGS1

while read BATCH_INFO_LINE
do
    read -r WGS_GROUP BATCH N_SAMPLES <<< "$BATCH_INFO_LINE"
    if [ $WGS_GROUP == "WGS7" ] || [ $WGS_GROUP == "WGS8" ] || [ $WGS_GROUP == "WGS9" ]; then continue; fi

    if [ $WGS_GROUP != $PREV_GROUP ]; then
        WGS_GROUP=$PREV_GROUP merge_WGS_GROUP
    fi

    PFX=$WGS_GROUP.batch$BATCH
    tar -xf /mnt/project/lohdata/david/mCAs_WGS/calls/$WGS_GROUP/$PFX.allelic_shift.tar -C ${TMP_DIR}

    echo chr{1..22} \
        | tr ' ' '\n' \
        | parallel --halt-on-error 2 --joblog=/dev/stderr WGS_GROUP=$WGS_GROUP PFX=$PFX extract_calls {}

    PREV_GROUP=$WGS_GROUP
done < /mnt/project/lohdata/resources/WGS_depth/sample_batches/WGS_group.batch.Nsamples.txt

WGS_GROUP=$PREV_GROUP merge_WGS_GROUP

echo chr{1..22} \
    | tr ' ' '\n' \
    | parallel --halt-on-error 2 --joblog=/dev/stderr "bcftools merge -m none $TMP_DIR/WGS{1..6}.{}.merged.bcf -Ob --no-index > $OUT_DIR/WGS_500k.{}.allelic_shift.bcf"
rm $TMP_DIR/WGS{1..6}.chr{1..22}.merged.bcf

echo chr{1..22} \
    | tr ' ' '\n' \
    | parallel --halt-on-error 2 --joblog=/dev/stderr "bcftools index -f $OUT_DIR/WGS_500k.{}.allelic_shift.bcf"

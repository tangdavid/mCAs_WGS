# set on input:
# CHR

set -euo pipefail

TMP_DIR=$HOME/tmp
mkdir -p $TMP_DIR

#source /mnt/project/lohdata/david/install_mocha.sh

ALLELIC_SHIFT_FILE=/mnt/project/lohdata/david/mCAs_WGS/calls/allelic_shift/WGS_500k.$CHR.allelic_shift.bcf
ls /mnt/project/lohdata/resources/burden_masks/vcfs/tmp/$CHR.LoF.scaffold.phased.*.bcf | sort -V > $TMP_DIR/concat_file.txt

bcftools concat -n -f $TMP_DIR/concat_file.txt -Ob -o $TMP_DIR/$CHR.LoF.scaffold.phased.bcf
bcftools index -f $TMP_DIR/$CHR.LoF.scaffold.phased.bcf

echo "`date`: annotating AS field into bcf..."
bcftools annotate -a $ALLELIC_SHIFT_FILE -c FMT/AS $TMP_DIR/$CHR.LoF.scaffold.phased.bcf -Ou -o $OUT_DIR/$CHR.LoF.scaffold.phased.AS.bcf
bcftools index -f $OUT_DIR/$CHR.LoF.scaffold.phased.AS.bcf

#echo "`date`: extending AS field to phased sites..."
#bcftools +extendFMT --format AS --phase --dist 500000 $TMP_DIR/$CHR.AS.bcf -Ob -o $OUT_BCF

#bcftools index -f $OUT_BCF

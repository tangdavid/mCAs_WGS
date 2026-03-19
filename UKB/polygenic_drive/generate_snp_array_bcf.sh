# set on input:
# OUT_DIR
set -euo pipefail

TMP_DIR=$HOME/tmp
mkdir -p $TMP_DIR

CHR=$1
bcftools view \
    -S <( bcftools query -l /mnt/project/lohdata/david/mCAs_WGS/calls/allelic_shift/WGS_500k.$CHR.allelic_shift.bcf) \
    /mnt/project/lohdata/david/mCAs_WGS/scaffold/$CHR.GRCh38.40709.AC.bcf \
    -Ob -o $TMP_DIR/$CHR.snp_array.calls_only.bcf 

bcftools index -f $TMP_DIR/$CHR.snp_array.calls_only.bcf

bcftools annotate \
    --columns FMT/AS \
    --annotations /mnt/project/lohdata/david/mCAs_WGS/calls/allelic_shift/WGS_500k.$CHR.allelic_shift.bcf \
    $TMP_DIR/$CHR.snp_array.calls_only.bcf \
    -Ob -o $OUT_DIR/$CHR.snp_array.AS.bcf 

bcftools index -f $OUT_DIR/$CHR.snp_array.AS.bcf

rm -f $TMP_DIR/$CHR.snp_array.calls_only.bcf{,.csi}

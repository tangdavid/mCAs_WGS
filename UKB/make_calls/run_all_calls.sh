# set on input
# WGS_GROUP
# BATCH
# OUT_DIR

set -euo pipefail

set +e 
sudo apt-get update
set -e
sudo apt --yes install parallel
sudo apt --yes install libdeflate-dev


PFX=${WGS_GROUP}.batch${BATCH}


TMP_DIR=$HOME/tmp
mkdir -p $TMP_DIR

echo "`date`: downloading CNV masks..."
zcat /mnt/project/lohdata/resources/WGS_depth/${WGS_GROUP}/${PFX}.CNVs.txt.gz \
    | grep -v MEDIAN \
    | cut -f1-5 \
    | cat - <(awk '{print "BLACKLIST\t"$0"\tMASK"}' /mnt/project/lohdata/resources/WGS_depth/lc_sv_cnv_masks.hg38.bed) \
    > ${TMP_DIR}/${PFX}.CNVs.txt

echo "`date`: downloading ref bias..."
cp /mnt/project/lohdata/david/mCAs_WGS/cramADs_WGS/10k.WGS.ref_bias.txt ${TMP_DIR}
chmod 644 ${TMP_DIR}/10k.WGS.ref_bias.txt

for CHR in chr{1..22}; do
    DX_PROJECT="WGS_500K"
    dx download "${DX_PROJECT}:/lohdata/david/mCAs_WGS/cramADs_WGS/${WGS_GROUP}/${PFX}.$CHR.correctedAD.bcf" -f -o $TMP_DIR
    bcftools index -f $TMP_DIR/${PFX}.$CHR.correctedAD.bcf
    ../../bin/mCAs_WGS make-calls \
        --cnv-mask ${TMP_DIR}/${PFX}.CNVs.txt \
        --ref-bias ${TMP_DIR}/10k.WGS.ref_bias.txt \
        --calls-out ${TMP_DIR}/${PFX}.$CHR.all_calls.txt \
        --bcf-out /dev/null \
        --bcf $TMP_DIR/${PFX}.$CHR.correctedAD.bcf
    rm $TMP_DIR/${PFX}.$CHR.correctedAD.bcf{,.csi}
done

awk 'NR==1 || FNR>1' ${TMP_DIR}/${PFX}.chr{1..22}.all_calls.txt > $OUT_DIR/$PFX.all_calls.txt 

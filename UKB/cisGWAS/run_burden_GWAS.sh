# set on input:
# CHR
# OUT_DIR
# REMOVE_LOCI
# CN

set -euo pipefail

set +e 
sudo apt-get update 2> /dev/null 1>&2
set -e
sudo apt --yes install libdeflate-dev 2> /dev/null 1>&2

TMP_DIR=$HOME/tmp
mkdir -p $TMP_DIR

cat /mnt/project/lohdata/resources/BOLT-LMM/w40709_CURRENT.FID_IID.txt \
    | awk '{print $1}' \
    > $TMP_DIR/remove.samples.txt 

CHR=$CHR bash remove_related_individuals.sh >> $TMP_DIR/remove.samples.txt

if [[ -v REMOVE_LOCI ]]; then 
    echo "`date`: removing CN-LOH explained by rare variants at ${REMOVE_LOCI} loci..."
    REMOVE_LOCI=$REMOVE_LOCI CHR=$CHR bash remove_rare_variant_carriers.sh >> $TMP_DIR/remove.samples.txt 
fi

echo "`date`: running GWAS for $CHR..."
INPUT_VCF=/mnt/project/lohdata/resources/burden_masks/vcfs/$CHR.LoF.scaffold.phased.AS.bcf
OUT_FILE=$OUT_DIR/$CHR.LoF.protein_coding.$CN.GWAS.txt.gz
CALLS_FILE=/mnt/project/lohdata/david/mCAs_WGS/calls/WGS_500k.calls.txt

../../bin/mCAs_WGS compute-cis-fisher \
    --calls $CALLS_FILE \
    --remove $TMP_DIR/remove.samples.txt \
    --cn $CN \
    < $INPUT_VCF | gzip > $OUT_FILE


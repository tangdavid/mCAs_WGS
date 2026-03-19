set -euo pipefail

CN=${CN:-"CN-LOH"}
DOWNSAMPLE=${DOWNSAMPLE:-1}

TMP_DIR=$HOME/tmp
CONDITIONAL=${CONDITIONAL:-1}
if [[ $CN != "CN-LOH" ]]; then CONDITIONAL=0; fi

mkdir -p $TMP_DIR

CHR=`echo $CHR_ARM | sed -E 's/chr([0-9,X]+)[pq]+/chr\1/g'`

if [[ $CHR != "chrX" ]]; then 
    CALLS_FILE=/mnt/project/lohdata/david/mCAs_WGS/calls/WGS_500k.calls.txt
    AS_BCF=/mnt/project/lohdata/david/mCAs_WGS/calls/allelic_shift/WGS_500k.$CHR.allelic_shift.bcf
else
    CALLS_FILE=/mnt/project/lohdata/david/mCAs_WGS/mLOX_new_phase/WGS_500k.mLOX.calls.txt 
    AS_BCF=/mnt/project/lohdata/david/mCAs_WGS/mLOX_new_phase/WGS_500k.mLOX.allelic_shift.bcf
fi

[ ! -f gcta64 ] && cp /mnt/project/lohdata/david/resources/gcta64 ./ && chmod 744 gcta64

[ ! -f $TMP_DIR/$CHR.prune.out ] && plink2 \
    --bcf $AS_BCF \
    --geno 0.05 \
    --indep-pairwise 500kb 0.5 \
    --out $TMP_DIR/$CHR \
    --set-all-var-ids @:#:\$r:\$a \
    --rm-dup exclude-mismatch \
    --new-id-max-allele-len 100 missing

cat /mnt/project/lohdata/resources/BOLT-LMM/{remove.nonEUR.FID_IID.40709.txt,w40709_CURRENT.FID_IID.txt} \
    | awk '{print $1}' \
    > $TMP_DIR/remove.samples.txt 

CN=$CN CHR=$CHR BCF_FILE=$AS_BCF CALL_FILE=$CALLS_FILE bash ../cisGWAS/remove_related_individuals.sh >> $TMP_DIR/remove.samples.txt

if [[ $CONDITIONAL == 1 ]]; then 
    echo "Removing CN-LOH on $CHR_ARM explained by rare variants..."
    CHR=$CHR bash ../cisGWAS/remove_rare_variant_carriers.sh >> $TMP_DIR/remove.samples.txt 
fi

bcftools query -l $AS_BCF \
    | awk -v mod=$DOWNSAMPLE 'NR%mod != 0' >> $TMP_DIR/remove.samples.txt

../../bin/write_CNLOH_GRM \
    --calls $CALLS_FILE \
    --chrom-arm $CHR_ARM \
    --grm-out $TMP_DIR/$CHR_ARM.gcta \
    --pheno-out $TMP_DIR/$CHR_ARM.gcta.pheno \
    --discard $TMP_DIR/$CHR.prune.out \
    --remove $TMP_DIR/remove.samples.txt \
    --cn $CN \
    < $AS_BCF

if [ ! -f $TMP_DIR/$CHR_ARM.gcta.pheno ]; then
    exit
fi

./gcta64 \
    --reml \
    --grm $TMP_DIR/$CHR_ARM.gcta \
    --pheno $TMP_DIR/$CHR_ARM.gcta.pheno \
    --reml-no-constrain \
    --prevalence 0.5 \
    --out $OUT_DIR/$CHR_ARM.gcta.$CN

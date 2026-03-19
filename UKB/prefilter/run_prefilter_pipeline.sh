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

source /mnt/project/lohdata/david/install_mocha.sh

PFX=${WGS_GROUP}.batch${BATCH}


TMP_DIR=$HOME/tmp
mkdir -p $TMP_DIR

echo "`date`: downloading imputed BCF..."
rm -f ${TMP_DIR}/${PFX}.imputedCommon.bcf{,.csi}
cp /mnt/project/lohdata/resources/imputedCommon/${WGS_GROUP}/${PFX}.imputedCommon.bcf{,.csi} ${TMP_DIR}
chmod u+w ${TMP_DIR}/${PFX}.imputedCommon.bcf{,.csi} 

ID_FILE="/mnt/project/lohdata/resources/WGS_depth/sample_batches/ID_40709.WGS_group.batch.txt"

echo "`date`: extracting DRAGEN VCF ADs..."
set +e # OK if some jobs fail
awk -v WGS_GROUP=$WGS_GROUP -v BATCH=$BATCH '$2==WGS_GROUP && $3==BATCH {print $1}' $ID_FILE \
    | parallel --joblog /dev/stderr --retries 4 "TMP_DIR=$TMP_DIR bash ../../tools/RAP/extract_DRAGEN_VCF.sh {} > ${TMP_DIR}/{}.hets.txt"
set -e

echo "`date`: downloading CNV mask..."
zcat /mnt/project/lohdata/resources/WGS_depth/${WGS_GROUP}/${PFX}.CNVs.txt.gz \
    | grep -v MEDIAN \
    | cut -f1-5 \
    | cat - <(awk '{print "BLACKLIST\t"$0"\tMASK"}' /mnt/project/lohdata/resources/WGS_depth/lc_sv_cnv_masks.hg38.bed) \
    > ${TMP_DIR}/${PFX}.CNVs.txt

echo "`date`: downloading ref bias file..."
cp /mnt/project/lohdata/david/mCAs_WGS/cramADs_WGS/10k.WGS.ref_bias.txt ${TMP_DIR}
chmod u+w ${TMP_DIR}/10k.WGS.ref_bias.txt

echo "`date`: creating input BCF..."
../../bin/mCAs_WGS prepare-bcf \
    --clear-ADs \
    --write-ADs \
    --het-format-str "${TMP_DIR}/%s.hets.txt" \
    --cnv-mask ${TMP_DIR}/${PFX}.CNVs.txt \
    --ref-bias-in ${TMP_DIR}/10k.WGS.ref_bias.txt \
    < ${TMP_DIR}/${PFX}.imputedCommon.bcf \
    > ${TMP_DIR}/${PFX}.VCF_AD.masked.bcf
bcftools index -f ${TMP_DIR}/${PFX}.VCF_AD.masked.bcf
rm -f ${TMP_DIR}/*.hets.txt

echo "`date`: running MoChA for prefiltering..."
bcftools +mocha \
    -g GRCh38 \
    -o ${TMP_DIR}/${PFX}.pre.bcf \
    -Ob \
    -c ${OUT_DIR}/${PFX}.pre.calls.tsv \
    -z ${OUT_DIR}/${PFX}.pre.stats.tsv \
    --LRR-GC-order 0 \
    --LRR-weight 0 \
    --bdev-LRR-BAF 6 \
    --auto-tel-pl 20 \
    --xy-major-pl 50 \
    --xy-minor-pl 40 \
    --flip-pl 40 \
    --mhc chr6:27518932-33480487 \
    --kir chr19:54071493-54992731 \
    --min-dist 50 \
    --threads 4 \
    ${TMP_DIR}/${PFX}.VCF_AD.masked.bcf

echo "`date`: downloading GRCh38 reference..."
rm -f ${TMP_DIR}/GRCh38_full_analysis_set_plus_decoy_hla.fa{,.fai}
cp /mnt/project/lohdata/resources/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa{,.fai} ${TMP_DIR}
chmod u+w ${TMP_DIR}/GRCh38_full_analysis_set_plus_decoy_hla.fa{,.fai}


cut -f1,3 ${OUT_DIR}/${PFX}.pre.calls.tsv \
    | tail +2 \
    | sort -u \
    | awk -v TMP_DIR=${TMP_DIR} '$2!="chrX" {print $0 > TMP_DIR"/"$2".candidates.txt"}'

echo "`date`: extracting cram ADs for candidate chromosomes..."
set +e
cat ${TMP_DIR}/*.candidates.txt | tr '\t' ' ' \
    | parallel --colsep ' ' --joblog /dev/stderr --retries 4 \
        "TMP_DIR=$TMP_DIR NGS=WGS bash ../../tools/RAP/extract_cram_AD_RAP.sh {1} {2} ${TMP_DIR}/${PFX}.imputedCommon.bcf | cut -f-6 > ${TMP_DIR}/{2}.ADs.{1}.txt"
set -e 

rm -f ${TMP_DIR}/tar.files.txt
for CHR in chr{1..22}; do 
    if [[ ! -f ${TMP_DIR}/${CHR}.candidates.txt ]]; then continue; fi
    echo "`date`: writing BCF for chromosome ${CHR}..."
    bcftools view -S <(cut -f1 ${TMP_DIR}/${CHR}.candidates.txt) -r ${CHR} ${TMP_DIR}/${PFX}.imputedCommon.bcf \
        | bcftools filter -i "INFO/AC!=INFO/AN && INFO/AC>0" -Ob \
        | ../../bin/mCAs_WGS prepare-bcf \
            --clear-ADs \
            --write-ADs \
            --het-format-str "${TMP_DIR}/${CHR}.ADs.%s.txt" \
            > ${TMP_DIR}/${PFX}.${CHR}.candidates.nomask.bcf 
    bcftools index -f  ${TMP_DIR}/${PFX}.${CHR}.candidates.nomask.bcf 

    echo ${TMP_DIR}/${PFX}.${CHR}.candidates.nomask.bcf >> ${TMP_DIR}/tar.files.txt 
    echo ${TMP_DIR}/${PFX}.${CHR}.candidates.nomask.bcf.csi >> ${TMP_DIR}/tar.files.txt 
done

echo "`date`: creating output tar..."
cat ${TMP_DIR}/tar.files.txt \
    | xargs -I{} basename {} \
    | xargs tar -cvf ${OUT_DIR}/${PFX}.candidates.nomask.tar -C ${TMP_DIR}

echo "done!"

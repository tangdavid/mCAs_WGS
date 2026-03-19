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
W_FILE="/mnt/project/lohdata/resources/BOLT-LMM/w40709_CURRENT.FID_IID.txt"
SEGDUPS="/mnt/project/lohdata/david/mCAs_WGS/GRCh38_supp_files/segdups.bed.gz"

echo "`date`: finding variants in segdups..."
bcftools annotate --no-version -Ou \
    -a $SEGDUPS \
    -c CHROM,FROM,TO,JK \
    -h <(echo '##INFO=<ID=JK,Number=1,Type=Float,Description="Jukes Cantor">' )  \
    ${TMP_DIR}/${PFX}.imputedCommon.bcf \
	| bcftools +fill-tags --no-version -Ou -t ^Y,MT,chrY,chrM -- -t ExcHet,F_MISSING \
	| bcftools view --no-version -Ou -G \
	| bcftools annotate --no-version -o $TMP_DIR/${PFX}.xcl.bcf -Ob --write-index \
		-i 'FILTER!="." && FILTER!="PASS" || INFO/JK<.02 || INFO/ExcHet<1e-6 || INFO/F_MISSING>1-.97' \
		-x ^INFO/JK,^INFO/ExcHet,^INFO/F_MISSING

echo "`date`: extracting DRAGEN VCF ADs..."
awk -v WGS_GROUP=$WGS_GROUP -v BATCH=$BATCH '
ARGIND==1 {removed[$1]=1} 
ARGIND==2 && $2==WGS_GROUP && $3==BATCH && !($1 in removed) {print $1}' $W_FILE $ID_FILE \
    | tee $TMP_DIR/kept_samples.txt \
    | parallel --joblog /dev/stderr --retries 4 --halt-on-error 2 "TMP_DIR=$TMP_DIR bash ../../tools/RAP/extract_DRAGEN_VCF.sh {} > ${TMP_DIR}/{}.hets.txt"

echo "`date`: downloading CNV mask..."
zcat /mnt/project/lohdata/resources/WGS_depth/${WGS_GROUP}/${PFX}.CNVs.txt.gz \
    | grep -v MEDIAN \
    | cut -f1-5 \
    | cat - <(awk '{print "BLACKLIST\t"$0"\tMASK"}' /mnt/project/lohdata/resources/WGS_depth/lc_sv_cnv_masks.hg38.bed) \
    > ${TMP_DIR}/${PFX}.CNVs.txt


echo "`date`: creating input BCF with CNVs masked..."
bcftools view -S $TMP_DIR/kept_samples.txt ${TMP_DIR}/${PFX}.imputedCommon.bcf -Ou \
    | ../../bin/mCAs_WGS prepare-bcf \
    --clear-ADs \
    --write-ADs \
    --het-format-str "${TMP_DIR}/%s.hets.txt" \
    --cnv-mask ${TMP_DIR}/${PFX}.CNVs.txt \
    > ${TMP_DIR}/${PFX}.VCF_AD.masked.bcf
bcftools index -f ${TMP_DIR}/${PFX}.VCF_AD.masked.bcf

echo "`date`: creating input BCF without CNVs masked..."
bcftools view -S $TMP_DIR/kept_samples.txt ${TMP_DIR}/${PFX}.imputedCommon.bcf -Ou \
    | bcftools +mochatools --no-version -Ou -- -t GC -f /mnt/project/lohdata/resources/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa \
    | ../../bin/mCAs_WGS prepare-bcf \
    --clear-ADs \
    --write-ADs \
    --het-format-str "${TMP_DIR}/%s.hets.txt" \
    > ${TMP_DIR}/${PFX}.VCF_AD.bcf
bcftools index -f ${TMP_DIR}/${PFX}.VCF_AD.bcf
rm -f ${TMP_DIR}/*.hets.txt

cat /mnt/project/lohdata/david/mCAs_WGS/sample_data/ID.genetic_sex.txt \
    | awk -v OFS='\t' '
BEGIN {print "sample_id","computed_gender"} 
$2==1 {print $1,"M"} 
$2==0 {print $1,"F"}'> ${TMP_DIR}/mocha.sex.txt

echo "`date`: running MoChA with TOPmed settings..."
bcftools +mocha \
    -g GRCh38 \
    -o /dev/null \
    -Ob \
    -c ${OUT_DIR}/${PFX}.mocha.TOPmed.calls.txt \
    -z /dev/null \
    --LRR-GC-order 0 \
    --LRR-weight 0 \
    --bdev-LRR-BAF 6 \
    --mhc chr6:27518932-33480487 \
    --kir chr19:54071493-54992731 \
    --min-dist 1000 \
    --threads `nproc` \
    --input-stats ${TMP_DIR}/mocha.sex.txt \
    ${TMP_DIR}/${PFX}.VCF_AD.masked.bcf

awk '
NR==1 {
    for(i=1; i<=NF; i++) 
        f[$i]=i; 
    print $0;
} 
{
    if($f["n_hets"] < 2000) next;
    if($f["chrom"] == "chrX") next;
    if($f["rel_cov"] > 2.9) next;
    if($f["rel_cov"]>2.5 && $f["bdev"] > 0.16) next;
    if($f["lod_baf_phase"] < 5) next; 
    if($f["bdev"] <= 0) next;
    print $0;
}' $OUT_DIR/${PFX}.mocha.TOPmed.calls.txt \
    > $OUT_DIR/${PFX}.mocha.TOPmed.filtered.calls.txt 

gzip $OUT_DIR/${PFX}.mocha.TOPmed.filtered.calls.txt
gzip $OUT_DIR/${PFX}.mocha.TOPmed.calls.txt

echo "`date`: running MoChA with default settings..."
bcftools +mocha \
    -g GRCh38 \
    -o /dev/null \
    -Ob \
    -c ${OUT_DIR}/${PFX}.mocha.default.calls.txt \
    -z ${TMP_DIR}/${PFX}.mocha.default.stats.txt \
    --mhc chr6:27518932-33480487 \
    --kir chr19:54071493-54992731 \
    --threads `nproc` \
    --cnp /mnt/project/lohdata/david/mCAs_WGS/GRCh38_supp_files/cnps.bed \
    --variants ^$TMP_DIR/${PFX}.xcl.bcf \
    --input-stats ${TMP_DIR}/mocha.sex.txt \
    ${TMP_DIR}/${PFX}.VCF_AD.bcf

awk -F "\t" '
ARGIND==1 && FNR==1 {for (i=1; i<=NF; i++) f[$i] = i}
ARGIND==1 && FNR>1 && ($(f["call_rate"])<.97 || $(f["baf_auto"])>.03) {xcl[$(f["sample_id"])]++}
ARGIND==2 && FNR==1 {for (i=1; i<=NF; i++) g[$i] = i; print}
ARGIND==2 && FNR>1 {
    gender=$(g["computed_gender"]); len=$(g["length"]); bdev=$(g["bdev"]);
    rel_cov=$(g["rel_cov"]); lod_baf_phase=$(g["lod_baf_phase"]); lod_baf_conc=$(g["lod_baf_conc"])
}
ARGIND==2 && FNR>1 && !($(g["sample_id"]) in xcl) && $(g["type"])!~"^CNP" &&
    ( $g["chrom"] != "chrX" ) &&
    ( $(g["bdev_se"])!="nan" || lod_baf_phase!="nan" && lod_baf_phase>10.0 ) &&
    ( rel_cov<2.1 || bdev<0.05 || len>5e5 && bdev<0.1 && rel_cov<2.5 || len>5e6 && bdev<0.15 )' \
    ${TMP_DIR}/${PFX}.mocha.default.stats.txt \
    ${OUT_DIR}/${PFX}.mocha.default.calls.txt \
    > $OUT_DIR/$PFX.mocha.default.filtered.calls.txt
gzip ${OUT_DIR}/${PFX}.mocha.default.calls.txt
gzip $OUT_DIR/$PFX.mocha.default.filtered.calls.txt

echo "done!"

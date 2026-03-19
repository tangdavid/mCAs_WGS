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

echo "`date`: downloading prefiltered BCFs..."
tar -xvf /mnt/project/lohdata/david/mCAs_WGS/prefilter/${WGS_GROUP}/${PFX}.candidates.nomask.tar -C ${TMP_DIR}

echo "`date`: downloading CNV masks..."
zcat /mnt/project/lohdata/resources/WGS_depth/${WGS_GROUP}/${PFX}.CNVs.txt.gz \
    | grep -v MEDIAN \
    | cut -f1-5 \
    | cat - <(awk '{print "BLACKLIST\t"$0"\tMASK"}' /mnt/project/lohdata/resources/WGS_depth/lc_sv_cnv_masks.hg38.bed) \
    > ${TMP_DIR}/${PFX}.CNVs.txt

echo "`date`: downloading ref-bias..."
cp /mnt/project/lohdata/david/mCAs_WGS/cramADs_WGS/10k.WGS.ref_bias.txt ${TMP_DIR}
chmod u+w ${TMP_DIR}/10k.WGS.ref_bias.txt

rm -f ${TMP_DIR}/tar.files.txt
rm -f ${TMP_DIR}/${PFX}.calls.no_depth.txt
for CHR in chr{1..22}; do
    if [[ ! -f ${TMP_DIR}/${PFX}.${CHR}.candidates.nomask.bcf ]]; then continue; fi
    echo "`date`: making calls for chromosome ${CHR}..."
    ../../bin/mCAs_WGS make-calls \
        --cnv-mask ${TMP_DIR}/${PFX}.CNVs.txt \
        --ref-bias ${TMP_DIR}/10k.WGS.ref_bias.txt \
        --calls-out ${TMP_DIR}/${PFX}.${CHR}.calls.txt \
        --bcf-out ${TMP_DIR}/${PFX}.${CHR}.allelic_shift.bcf \
        --bcf ${TMP_DIR}/${PFX}.${CHR}.candidates.nomask.bcf
    bcftools index -f ${TMP_DIR}/${PFX}.${CHR}.allelic_shift.bcf

    if [[ ! -f ${TMP_DIR}/${PFX}.calls.no_depth.txt ]]; then 
        cat ${TMP_DIR}/${PFX}.${CHR}.calls.txt > ${TMP_DIR}/${PFX}.calls.no_depth.txt
    else 
        awk 'FNR!=1' ${TMP_DIR}/${PFX}.${CHR}.calls.txt >> ${TMP_DIR}/${PFX}.calls.no_depth.txt
    fi
    echo ${TMP_DIR}/${PFX}.${CHR}.allelic_shift.bcf >> ${TMP_DIR}/tar.files.txt
    echo ${TMP_DIR}/${PFX}.${CHR}.allelic_shift.bcf.csi >> ${TMP_DIR}/tar.files.txt
done 


echo "`date`: annotating in PCadj depths and sex..."

awk -v OFS='\t' 'FNR!=1 {print $1,$3":"$4"-"$5}' ${TMP_DIR}/${PFX}.calls.no_depth.txt > ${TMP_DIR}/${PFX}.ID_region_list.txt
ID_REGION_LIST=${TMP_DIR}/${PFX}.ID_region_list.txt OUT_FILE=${TMP_DIR}/${PFX}.depth.txt bash ../../tools/RAP/compute_region_depth_PCadj_RAP.sh

SEX_FILE="/mnt/project/lohdata/resources/BOLT-LMM/ID_40709.cov_SEX_GENETIC.txt"

awk -v OFS='\t' '
    ARGIND==1 && FNR==1 { for (i=1; i<=NF; i++) f1[$i]=i; }
    ARGIND==1 && FNR>1 {
        ID=$(f1["ID"]);
        chrStr=$(f1["chr"]);
        bpStart=$(f1["bpStart"]);
        bpEnd=$(f1["bpEnd"]);

        depth[ID,chrStr,bpStart,bpEnd]=$5; 
        depth_se[ID,chrStr,bpStart,bpEnd]=$6
    } 
    ARGIND==2 { 
        sex[$1] = ($2==1) ? "M" : "F";
    }
    ARGIND==3 && FNR==1 { 
        for (i=1; i<=NF; i++) f3[$i]=i; 
        print $0;
    }
    ARGIND==3 && FNR>1 {
        ID=$(f3["ID"]);
        chrStr=$(f3["chr"]);
        bpStart=$(f3["bpStart"]);
        bpEnd=$(f3["bpEnd"]);

        $(f3["sex"])=sex[ID]; 
        $(f3["depth"])=depth[ID,chrStr,bpStart,bpEnd]; 
        $(f3["depthSE"])=depth_se[ID,chrStr,bpStart,bpEnd]; 
        print $0
    }' ${TMP_DIR}/${PFX}.depth.txt ${SEX_FILE} ${TMP_DIR}/${PFX}.calls.no_depth.txt  \
    > ${TMP_DIR}/${PFX}.calls.without_isochromosomes.txt 

# call isochromosome events
WGS_GROUP=$WGS_GROUP \
BATCH=$BATCH \
ID_REGION_LIST=<(awk 'NR>1 {print $1,$3":"1"-"1000000000}' ${TMP_DIR}/${PFX}.calls.without_isochromosomes.txt) \
OUT_FILE=/dev/null \
    bash ../../tools/RAP/compute_region_depth_PCadj_RAP.sh /dev/stdout \
    | BIN_SIZE=10 bash ../../tools/bin_depth_profile.sh \
    | awk -v OFS='\t' 'NR==1{$2="chr"} NR>1 {$2=substr($2,1,index($2,":")-1)} {print}' \
    > $TMP_DIR/$PFX.depth_profiles.txt 

CALLS_FILE=${TMP_DIR}/${PFX}.calls.without_isochromosomes.txt \
DEPTH_PROFILES=${TMP_DIR}/$PFX.depth_profiles.txt \
    bash call_isochromosomes.sh | gzip > ${OUT_DIR}/${PFX}.calls.txt.gz

for CHR in chr{1..22}; do
    echo "`date`: validating calls on ${CHR}..."
    TMP_DIR=${TMP_DIR} bash validate_calls_WES.sh \
        ${TMP_DIR}/${PFX}.${CHR}.calls.txt \
        ${TMP_DIR}/${PFX}.${CHR}.allelic_shift.bcf \
        ${TMP_DIR}/${PFX}.${CHR}.validation.WES.txt 
done

echo -e "ID\tchr\tbpStart\tbpEnd\toverRepWES\tunderRepWES\tzscore\tbdev\tbdevSE\texpectedValRate" \
    > ${OUT_DIR}/${PFX}.validation.WES.txt
cat ${TMP_DIR}/${PFX}.chr{1..22}.validation.WES.txt >> ${OUT_DIR}/${PFX}.validation.WES.txt
gzip ${OUT_DIR}/${PFX}.validation.WES.txt

echo "`date`: creating output tar..."
cat ${TMP_DIR}/tar.files.txt \
    | xargs -I{} basename {} \
    | xargs tar -cvf ${OUT_DIR}/${PFX}.allelic_shift.tar -C ${TMP_DIR}

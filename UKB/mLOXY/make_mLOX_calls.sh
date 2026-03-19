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

#CRAM_AD_BCF="/mnt/project/lohdata/david/mCAs_WGS/chrX_ADs/${WGS_GROUP}/${PFX}.chrX.correctedAD.bcf"
CRAM_AD_BCF="/mnt/project/lohdata/yichen/data/mloxy/mlox/phased_ad/${WGS_GROUP}/${PFX}.chrX.correctedAD.newphase.bcf"

TMP_DIR=$HOME/tmp
mkdir -p $TMP_DIR

zcat /mnt/project/lohdata/resources/WGS_depth/${WGS_GROUP}/${PFX}.CNVs.txt.gz \
    | grep -v MEDIAN \
    | cut -f1-5 \
    | cat - <(awk '{print "BLACKLIST\t"$0"\tMASK"}' /mnt/project/lohdata/resources/WGS_depth/lc_sv_cnv_masks.hg38.bed ) \
    > ${TMP_DIR}/${PFX}.CNVs.txt

cp /mnt/project/lohdata/david/mCAs_WGS/cramADs_WGS/10k.WGS.ref_bias.txt ${TMP_DIR}
chmod 644 ${TMP_DIR}/10k.WGS.ref_bias.txt

../../bin/mCAs_WGS make-calls \
    --cnv-mask ${TMP_DIR}/${PFX}.CNVs.txt \
    --ref-bias ${TMP_DIR}/10k.WGS.ref_bias.txt \
    --calls-out ${TMP_DIR}/${PFX}.mLOX.calls.no_depth.txt \
    --bcf-out ${OUT_DIR}/${PFX}.mLOX.allelic_shift.bcf \
    --bcf $CRAM_AD_BCF
bcftools index -f ${OUT_DIR}/${PFX}.mLOX.allelic_shift.bcf

sed -i 's/chr23/chrX/g' ${TMP_DIR}/${PFX}.mLOX.calls.no_depth.txt

echo "`date`: annotating in PCadj depths and sex..."

awk -v OFS='\t' 'FNR!=1 {print $1,$3":"$4"-"$5}' ${TMP_DIR}/${PFX}.mLOX.calls.no_depth.txt > ${TMP_DIR}/${PFX}.ID_region_list.txt
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
    }' ${TMP_DIR}/${PFX}.depth.txt ${SEX_FILE} ${TMP_DIR}/${PFX}.mLOX.calls.no_depth.txt  \
    > ${OUT_DIR}/${PFX}.mLOX.calls.txt 

set +e
TMP_DIR=${TMP_DIR} bash ../make_calls/validate_calls_WES.sh \
    ${OUT_DIR}/${PFX}.mLOX.calls.txt \
    ${OUT_DIR}/${PFX}.mLOX.allelic_shift.bcf \
    ${OUT_DIR}/${PFX}.mLOX.validation.WES.txt 
set -e 

gzip ${OUT_DIR}/${PFX}.mLOX.validation.WES.txt
gzip ${OUT_DIR}/${PFX}.mLOX.calls.txt

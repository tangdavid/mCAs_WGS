# set on input:
# FIRST_BATCH
# OUT_DIR
# FOCAL_REGION

set -euo pipefail

set +e 
sudo apt-get update
set -e
sudo apt --yes install parallel
sudo apt --yes install libdeflate-dev
sudo apt --yes install pigz

REGIONS=`grep -w $FOCAL_REGION ./focal_region_list.txt | awk '{print $2}'`
echo $FOCAL_REGION $REGIONS

awk 'substr($1,4) <= 6' /mnt/project/lohdata/resources/WGS_depth/sample_batches/WGS_group.batch.Nsamples.txt \
    | awk -v FIRST_BATCH=$FIRST_BATCH 'FNR>=FIRST_BATCH*10 && FNR<FIRST_BATCH*10+10' \
    | while read LINE;  do 

    read -r WGS_GROUP BATCH N_SAMP <<< $LINE
    echo $WGS_GROUP $BATCH  >> /dev/stderr
    WGS_GROUP=$WGS_GROUP BATCH=$BATCH bash extract_focal_region_depth.sh $REGIONS \
        | awk -v OFS='\t' '
        ARGIND==1 {region2name[$2]=$1} 
        ARGIND==2 && FNR==1 {for (i=1; i<=NF; i++) f[$i]=i; print $0,"region"} 
        ARGIND==2 && FNR>1 {region=$f["chr"]":"$f["bpStart"]"-"$f["bpEnd"]; print $0,region2name[region]}' focal_region_list.txt - \
        | gzip > ${OUT_DIR}/$WGS_GROUP.batch$BATCH.$FOCAL_REGION.depth.txt.gz 
done 

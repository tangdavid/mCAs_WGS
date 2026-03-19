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

source /mnt/project/lohdata/david/tools/python_setup_RAP.sh

TMP_DIR=$HOME/tmp
mkdir -p $TMP_DIR

REGIONS=`grep -w $FOCAL_REGION ./focal_region_list.txt | awk '{print $2}'`

echo $FOCAL_REGION $REGIONS
awk 'substr($1,4) <= 6' /mnt/project/lohdata/resources/WGS_depth/sample_batches/WGS_group.batch.Nsamples.txt \
    | awk -v FIRST_BATCH=$FIRST_BATCH 'FNR>=FIRST_BATCH*10 && FNR<FIRST_BATCH*10+10' \
    | while read LINE;  do 

    read -r WGS_GROUP BATCH N_SAMP <<< $LINE
    echo $WGS_GROUP $BATCH  >> /dev/stderr
    WGS_GROUP=$WGS_GROUP BATCH=$BATCH bash find_focal_breakpoint.sh $REGIONS \
        | sort -k1,1n -k2,2V -k3,3n -k4,4n -u \
        | gzip > ${OUT_DIR}/$WGS_GROUP.batch$BATCH.$FOCAL_REGION.depth_breakpoint.txt.gz 
done 

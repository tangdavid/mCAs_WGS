# set on input:
# WGS_GROUP
# BATCH

set -euo pipefail

TMP_DIR=$HOME/tmp
mkdir -p $TMP_DIR

REGIONS=`echo ${@} | tr ' ' ','`

awk -v WGS_GROUP=$WGS_GROUP -v BATCH=$BATCH -v REGIONS=$REGIONS '
BEGIN{split(REGIONS,regions,",")}

$2==WGS_GROUP && $3==BATCH {
    split(REGIONS, regions, ",")
    for (i in regions) {
	    print $1,regions[i];	
    }
}
' /mnt/project/lohdata/resources/WGS_depth/sample_batches/ID_40709.WGS_group.batch.txt \
    > $TMP_DIR/$WGS_GROUP.batch$BATCH.ID_region.txt 


WGS_GROUP=$WGS_GROUP \
BATCH=$BATCH \
ID_REGION_LIST=$TMP_DIR/$WGS_GROUP.batch$BATCH.ID_region.txt \
OUT_FILE=/dev/stdout \
    bash ../../tools/RAP/compute_region_depth_PCadj_RAP.sh /dev/null

rm $TMP_DIR/$WGS_GROUP.batch$BATCH.ID_region.txt

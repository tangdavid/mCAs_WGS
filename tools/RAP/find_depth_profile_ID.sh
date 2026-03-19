# set on input:

set -euo pipefail

ID=$1
REGIONS=${@:2}
echo $REGIONS
read -r ID WGS_GROUP BATCH <<< `grep $ID /mnt/project/lohdata/resources/WGS_depth/sample_batches/ID_40709.WGS_group.batch.txt`

WGS_GROUP=$WGS_GROUP \
BATCH=$BATCH \
ID_REGION_LIST=<(echo $REGIONS | tr ' ' '\n' | awk -v ID=$ID -v OFS='\t' '{print ID, $1}' ) \
OUT_FILE=/dev/null \
    bash `dirname ${BASH_SOURCE[0]}`/compute_region_depth_PCadj_RAP.sh /dev/stdout 

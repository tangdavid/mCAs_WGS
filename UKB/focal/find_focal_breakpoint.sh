# set on input:
# WGS_GROUP
# BATCH

set -euo pipefail

MAX_LENGTH="${MAX_LENGTH:=5000}"
MIN_LENGTH="${MIN_LENGTH:=100}"

echo -e "MAX_LENGTH: $MAX_LENGTH\nMIN_LENGTH: $MIN_LENGTH" > /dev/stderr

REGIONS=${@}

TMP_DIR=$HOME/tmp
mkdir -p $TMP_DIR

rm -f $TMP_DIR/$WGS_GROUP.batch$BATCH.discordant.binned.txt

echo "downloading discordant reads..." > /dev/stderr
DISCORDANT_FILE=$TMP_DIR/$WGS_GROUP.batch$BATCH.discordant.txt
unpigz -c /mnt/project/lohdata/resources/WGS_depth/$WGS_GROUP/$WGS_GROUP.batch$BATCH.discordant.txt.gz > $DISCORDANT_FILE

for REGION in $REGIONS; do
    echo "binning discordant reads in $REGION..." > /dev/stderr
    REGION=$REGION \
    MAX_LENGTH=$MAX_LENGTH \
    MIN_LENGTH=$MIN_LENGTH \
    DISCORDANT_FILE=$DISCORDANT_FILE \
        bash bin_discordant.sh >> $TMP_DIR/$WGS_GROUP.batch$BATCH.discordant.binned.txt
done

if [ `cat $TMP_DIR/$WGS_GROUP.batch$BATCH.discordant.binned.txt | wc -l` == 0 ]; then
    echo "Warning: no regions in $TMP_DIR/$WGS_GROUP.batch$BATCH.discordant.binned.txt"  >> /dev/stderr 
    exit 0
fi

echo "extracting depth for discordant reads..." > /dev/stderr
WGS_GROUP=$WGS_GROUP \
BATCH=$BATCH \
ID_REGION_LIST=<(awk '{print $1,$5}' $TMP_DIR/$WGS_GROUP.batch$BATCH.discordant.binned.txt | awk '{split($2,a,":"); split(a[2],pos,"-"); window=a[1]":"pos[1]-5e6"-"pos[2]+5e6; print $1,window}' | sort -u) \
OUT_FILE=/dev/null \
    bash ../../tools/RAP/compute_region_depth_PCadj_RAP.sh /dev/stdout \
    | bash get_discordant_read_depth.sh $TMP_DIR/$WGS_GROUP.batch$BATCH.discordant.binned.txt -

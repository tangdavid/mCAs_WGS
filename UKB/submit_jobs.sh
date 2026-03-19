set -euo pipefail

LAUNCH_SCRIPT=$1
NAME=$2

while read BATCH_INFO_LINE
do
    read -r WGS_GROUP BATCH N_SAMPLES <<< "$BATCH_INFO_LINE"
    #if [ $WGS_GROUP == "WGS4" ]; then continue; fi
    if [ $WGS_GROUP == "WGS7" ]; then continue; fi
    if [ $WGS_GROUP == "WGS8" ]; then continue; fi
    if [ $WGS_GROUP == "WGS9" ]; then continue; fi
    WGS_GROUP=$WGS_GROUP BATCH=$BATCH NAME=$NAME bash $LAUNCH_SCRIPT 

done < ./sample_batches/WGS_group.batch.Nsamples.txt

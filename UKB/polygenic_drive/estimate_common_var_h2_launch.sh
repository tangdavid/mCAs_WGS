# set on input:
# CHR_ARM
# CN

set -euo pipefail

DIR="/lohdata/david/mCAs_WGS"
CN=${CN:-"CN-LOH"}
DOWNSAMPLE=${DOWNSAMPLE:-1}

if [[ $CHR_ARM == "chrXpq" ]]; then 
    INSTANCE_TYPE=mem3_ssd1_v2_x4
else
    INSTANCE_TYPE=mem1_ssd1_v2_x4
fi

echo "estimating common var h2 for $CHR_ARM $CN"
dx run swiss-army-knife \
    -iin="$DIR"/run_pipeline.sh \
    -icmd="set -euo pipefail && CHR_ARM=$CHR_ARM CN=$CN DOWNSAMPLE=$DOWNSAMPLE /usr/bin/time -v bash run_pipeline.sh polygenic_drive estimate_common_var_h2.sh 2>&1 | tee $CHR_ARM.$CN.h2.dx_run.log" \
    --instance-type $INSTANCE_TYPE \
    --destination /lohdata/david/mCAs_WGS/polygenic_drive/h2/ \
    --priority low \
    --brief \
    --ignore-reuse \
    --name estimate_common_var_h2.$CN.$CHR_ARM \
    --allow-ssh \
    --tag estimate_common_var_h2.$CN \
    -y

# set on input:
# CHR
# CONDITIONAL
# CN

set -euo pipefail

PRIORITY=${PRIORITY:-low}
echo -n "$PRIORITY priority; "

DIR="/lohdata/david/mCAs_WGS"

dx ls "/lohdata/resources/burden_masks/vcfs" | grep -q $CHR.LoF.scaffold.phased.AS.bcf || exit

if [[ $CONDITIONAL == 1 ]]; then
    DESTINATION="$DIR"/cisGWAS/results/conditional/
    echo "$CHR $CN conditional analysis"
else 
    DESTINATION="$DIR"/cisGWAS/results/
    echo "$CHR $CN"
fi

dx run /lohdata/resources/hts-knife \
    -iin="$DIR"/run_pipeline.sh \
    -icmd="set -euo pipefail && CN=$CN CONDITIONAL=$CONDITIONAL CHR=$CHR /usr/bin/time -v bash run_pipeline.sh cisGWAS run_burden_GWAS.sh 2>&1 | tee $CHR.$CN.burden_GWAS.dx_run.log" \
    --instance-type mem1_ssd1_v2_x2 \
    --destination $DESTINATION \
    --priority $PRIORITY \
    --brief \
    --ignore-reuse \
    --name $CHR.$CN.burden_GWAS \
    --allow-ssh \
    --tag $CN.burden_GWAS \
    -y

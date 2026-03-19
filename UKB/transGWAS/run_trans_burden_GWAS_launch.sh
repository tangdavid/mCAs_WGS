# set on input:

set -euo pipefail

DIR="/lohdata/david/mCAs_WGS"
PRIORITY=${PRIORITY:-low}
echo "Submitting with $PRIORITY priority"

for CHR in {1..22}; do
    echo chr$CHR
    dx run swiss-army-knife \
        -iin="$DIR"/run_pipeline.sh \
        -icmd="set -euo pipefail && CHR=$CHR /usr/bin/time -v bash run_pipeline.sh transGWAS run_trans_burden_GWAS.sh 2>&1 | tee burden.transGWAS.c$CHR.dx_run.log" \
        --instance-type mem1_ssd1_v2_x4 \
        --destination "$DIR"/transGWAS/burden_results \
        --priority $PRIORITY \
        --brief \
        --ignore-reuse \
        --name burden.transGWAS.c$CHR \
        --tag burden.transGWAS \
        --allow-ssh \
        -y
done

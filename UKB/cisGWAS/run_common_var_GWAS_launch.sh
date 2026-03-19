# set on input:
# CHR
# REMOVE_LOCI

set -euo pipefail

DIR="/lohdata/david/mCAs_WGS"

if [[ $REMOVE_LOCI == "sig" ]]; then 
    DESTINATION="$DIR"/cisGWAS/results/common_var/restricted/
    echo $CHR without rare variants
    ICMD="set -euo pipefail && REMOVE_LOCI=$REMOVE_LOCI CHR=$CHR /usr/bin/time -v bash run_pipeline.sh cisGWAS run_common_var_GWAS.sh 2>&1 | tee $CHR.common_var_GWAS.dx_run.log"
else 
    DESTINATION="$DIR"/cisGWAS/results/common_var/
    echo $CHR
    ICMD="set -euo pipefail && CHR=$CHR /usr/bin/time -v bash run_pipeline.sh cisGWAS run_common_var_GWAS.sh 2>&1 | tee $CHR.common_var_GWAS.dx_run.log"
fi

dx run /lohdata/resources/hts-knife \
    -iin="$DIR"/run_pipeline.sh \
    -icmd="$ICMD" \
    --instance-type mem1_ssd1_v2_x4 \
    --destination "$DESTINATION" \
    --priority low \
    --brief \
    --ignore-reuse \
    --name $CHR.common_var_GWAS \
    --allow-ssh \
    --tag common_var_GWAS \
    -y

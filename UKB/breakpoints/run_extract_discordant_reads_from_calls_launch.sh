# set on input:

set -euo pipefail

DIR="/lohdata/david/mCAs_WGS"

CORES=8

dx run /lohdata/resources/hts-knife \
    -iin="$DIR"/run_pipeline.sh \
    -icmd="set -euo pipefail && /usr/bin/time -v bash run_pipeline.sh breakpoints run_extract_discordant_reads_from_calls.sh 2>&1 | tee extract_discordant_reads.dx_run.log" \
    --instance-type mem1_ssd1_v2_x$CORES \
    --destination "$DIR"/breakpoints/ \
    --priority high \
    --brief \
    --ignore-reuse \
    --name extract_discordant_reads \
    --allow-ssh \
    -y

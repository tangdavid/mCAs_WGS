set -euo pipefail

DIR="/lohdata/david/mCAs_WGS"

CORES=4

for CHR in {1..22} X
do
    dx run /lohdata/resources/hts-knife \
	-iin="${DIR}/run_pipeline.sh" \
	-icmd="set -euo pipefail && CHR=$CHR /usr/bin/time -v bash run_pipeline.sh generate_impute_ref run_generate_reference_RAP.sh 2>&1 | tee $CHR.dx_run.log" \
	--instance-type mem2_ssd1_v2_x$CORES \
	--destination "${DIR}/imputation_reference/" \
	--priority low \
	--brief \
	--ignore-reuse \
	--name generate_imputation_reference.chr$CHR \
	-y
done

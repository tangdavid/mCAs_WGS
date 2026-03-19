set -euo pipefail

DIR="/lohdata/david/mCAs_WGS"

CORES=2

for CHR in {1..22} X
do
    dx run /lohdata/resources/hts-knife \
    -iin="/lohdata/resources/mosaic500K_phasing/chr${CHR}.GRCh38.40709.bcf" \
	-icmd="set -euo pipefail && bcftools view chr${CHR}.GRCh38.40709.bcf | bcftools +fill-tags -- -t AC,AN | bcftools view -Ob -o chr${CHR}.GRCh38.40709.AC.bcf && bcftools index -f chr${CHR}.GRCh38.40709.AC.bcf 2>&1 | tee $CHR.dx_run.log" \
	--instance-type mem1_ssd1_v2_x$CORES \
	--destination "${DIR}/scaffold/" \
	--priority low \
	--brief \
	--ignore-reuse \
	--name add_AC_scaffold.chr$CHR \
	-y
done

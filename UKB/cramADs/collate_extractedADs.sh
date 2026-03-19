# set on input:
# WGS_GROUP
# BATCH
# NGS
# TMP_DIR

set -euo pipefail

PFX=${WGS_GROUP}.batch${BATCH}

bcftools concat \
    -n /mnt/project/lohdata/david/mCAs_WGS/cramADs_${NGS}/${WGS_GROUP}/${PFX}.chr{{1..22},X}.correctedAD.bcf \
    -o ${TMP_DIR}/${PFX}.correctedAD.bcf
bcftools index -f ${TMP_DIR}/${PFX}.correctedAD.bcf

bcftools concat \
    -n /mnt/project/lohdata/david/mCAs_WGS/cramADs_${NGS}/${WGS_GROUP}/${PFX}.chr{{1..22},X}.rawAD.bcf \
    -o ${TMP_DIR}/${PFX}.rawAD.bcf
bcftools index -f ${TMP_DIR}/${PFX}.rawAD.bcf

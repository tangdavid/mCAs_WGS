# set on input: $CHR

set -euo pipefail

set +e # OK if apt-get update fails
sudo apt-get update
set -e
sudo apt --yes install libdeflate-dev 

unset DX_WORKSPACE_ID # to allow dx make_download_url

TMP_DIR=$HOME/tmp
mkdir -p $TMP_DIR


ID_FILE="/mnt/project/lohdata/ploh/mCAs_WGS/sample_batches/ID_40709.WGS_group.batch.txt"
PHASED_BCF="/mnt/project/Bulk/Whole genome sequences/SHAPEIT Phased VCFs/ukb20279_c${CHR}_b0_v1.vcf.gz"

# create sample lists for both references

awk '$2=="WGS3" {print $1}' $ID_FILE \
    > $TMP_DIR/WGS3_IDs.txt

awk '$2!="WGS3" {print $1}' $ID_FILE \
    > $TMP_DIR/WGS1_WGS2_IDs.txt

../../bin/vcf_extract_common_AC_AN "${PHASED_BCF}" \
    | bcftools view -Ou -o $TMP_DIR/chr${CHR}.phased.bcf 

echo "done extracting common variants"


for REF in WGS1_WGS2 WGS3; do
    bcftools view $TMP_DIR/chr${CHR}.phased.bcf \
        -S $TMP_DIR/${REF}_IDs.txt \
        -Ob -o $TMP_DIR/chr${CHR}.${REF}.bcf \
        && bcftools index -f $TMP_DIR/chr${CHR}.${REF}.bcf 

    mkdir ${REF}
    
    ../../bin/xcftools_static view \
        -i $TMP_DIR/chr${CHR}.${REF}.bcf \
        --o $REF/chr${CHR}.phased.${REF}.xcf.bcf \
        -O sh \
        --region chr${CHR} \
        -T 4
    
    rm $TMP_DIR/chr${CHR}.${REF}.bcf*
    echo "done creating ref panel ${REF}"
done

echo "done!"

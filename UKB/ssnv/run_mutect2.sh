set -euo pipefail

CRAM_DIR="WGS_500K:/Bulk/DRAGEN WGS/Whole genome CRAM files (DRAGEN) [500k release]/${ID:0:2}"

REF_FILE=$TMP_DIR/GRCh38_full_analysis_set_plus_decoy_hla.fa

REGION=chr13:48000000-52000000

GATK_DIR=$TMP_DIR/gatk-4.3.0.0
SNPEFF_DIR=$TMP_DIR/snpEff
GNOMAD_DIR=$TMP_DIR/gnomad

$GATK_DIR/gatk Mutect2 \
    -R $REF_FILE \
    -I `dx make_download_url "$CRAM_DIR/${ID}_24048_0_0.dragen.cram"` \
    --read-index `dx make_download_url "$CRAM_DIR/${ID}_24048_0_0.dragen.cram.crai"` \
    -O $TMP_DIR/$ID.vcf.gz \
    --germline-resource $GNOMAD_DIR/af-only-gnomad.hg38.vcf.gz \
    -L $REGION

$GATK_DIR/gatk FilterMutectCalls \
    -V $TMP_DIR/$ID.vcf.gz \
    -R $REF_FILE \
    -O $TMP_DIR/$ID.filtered.vcf.gz

rm $TMP_DIR/$ID.vcf.gz*

bcftools view $TMP_DIR/$ID.filtered.vcf.gz -i 'FILTER=="PASS"' \
    | java -Xmx4g -jar $SNPEFF_DIR/snpEff.jar GRCh38.p7.RefSeq \
    | bcftools view -i 'INFO/ANN ~ "HIGH"' \
    > $OUT_DIR/$ID.high_vep.vcf

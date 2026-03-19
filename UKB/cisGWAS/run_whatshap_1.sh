# set on input:
# ID
# CHR
# POS
# REFERENCE

set -euo pipefail 

REGION=$CHR:$[POS-1000]-$[POS+1000]

TMP_DIR=$HOME/tmp
mkdir -p $TMP_DIR


CRAM_DIR="/Bulk/DRAGEN WGS/Whole genome CRAM files (DRAGEN) [500k release]/${ID:0:2}"
CRAM_URL=`dx make_download_url "${CRAM_DIR}/${ID}_24048_0_0.dragen.cram"`
CRAI_URL=`dx make_download_url "${CRAM_DIR}/${ID}_24048_0_0.dragen.cram.crai"`

#VCF_FILE="/mnt/project/Bulk/DRAGEN WGS/Whole genome variant call files (VCFs) (DRAGEN) [500k release]/${ID:0:2}/${ID}_24053_0_0.dragen.hard-filtered.vcf.gz"
VCF_DIR="/Bulk/DRAGEN WGS/Whole genome variant call files (VCFs) (DRAGEN) [500k release]/${ID:0:2}"
VCF_URL=`dx make_download_url "${VCF_DIR}/${ID}_24053_0_0.dragen.hard-filtered.vcf.gz"`
TBI_URL=`dx make_download_url "${VCF_DIR}/${ID}_24053_0_0.dragen.hard-filtered.vcf.gz.tbi"`

PREFIX="$ID."`echo $REGION | tr ':-' '_'`

bcftools view -r $REGION "$VCF_URL"##idx##$TBI_URL -g het -m2 -M2 > $TMP_DIR/$PREFIX.input.vcf
rm -f ${ID}_24053_0_0.dragen.hard-filtered.vcf.gz.tbi
samtools view -T $REFERENCE $CRAM_URL -X $CRAI_URL -b $REGION > $TMP_DIR/$PREFIX.input.bam
samtools index $TMP_DIR/$PREFIX.input.bam 

whatshap phase -o /dev/stdout $TMP_DIR/$PREFIX.input.vcf $TMP_DIR/$PREFIX.input.bam --reference $REFERENCE \
    | bcftools query -f '[%CHROM:%POS:%REF:%ALT\t%GT\t%PS\n]' \
    | awk -v OFS='\t' -v ID=$ID '{print ID,$0}'

rm $TMP_DIR/$PREFIX.input.{bam,bam.bai,vcf}


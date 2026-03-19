# set on input:
# ID
# POS
# CHR
# REFERENCE

set -euo pipefail

TMP_DIR=$HOME/tmp
mkdir -p $TMP_DIR

REGION=$CHR:$[POS-1000]-$[POS+1000]

DX_PROJECT="WGS_500K"
CRAM_DIR="${DX_PROJECT}:/Bulk/DRAGEN WGS/Whole genome CRAM files (DRAGEN) [500k release]/${ID:0:2}"
CRAM_URL=`dx make_download_url "${CRAM_DIR}/${ID}_24048_0_0.dragen.cram"`

VCF_DIR="${DX_PROJECT}:/Bulk/DRAGEN WGS/Whole genome variant call files (VCFs) (DRAGEN) [500k release]/${ID:0:2}"
VCF_URL=`dx make_download_url "${VCF_DIR}/${ID}_24053_0_0.dragen.hard-filtered.vcf.gz"`
TBI_URL=`dx make_download_url "${VCF_DIR}/${ID}_24053_0_0.dragen.hard-filtered.vcf.gz.tbi"`

PREFIX="$ID."`echo $REGION | tr ':-' '_'`
bcftools view -r $REGION "$VCF_URL"##idx##$TBI_URL -g het -m2 -M2 > $TMP_DIR/$PREFIX.input.vcf
rm -f ${ID}_24053_0_0.dragen.hard-filtered.vcf.gz.tbi
CRAI_URL=`dx make_download_url "${CRAM_DIR}/${ID}_24048_0_0.dragen.cram.crai"`

samtools mpileup \
    -Q 20 \
    --output-QNAME \
    -l $TMP_DIR/$PREFIX.input.vcf \
    -f $REFERENCE \
    -r $CHR:$((POS-1000))-$((POS+1000)) \
    $CRAM_URL -X $CRAI_URL \
    | cut -f2,4,5,7 \
    | sed 's/\^.//g' \
    | tr -d '$' \
    | perl -lane '$F[2] =~ s/(.([+-](\d+)(??{".{$3}"}))*)/ $1/g; print "@F"' \
    | awk -v OFS='\t' -v POS=$POS -v ID=$ID '
ARGIND==1 {
    variant2pos[$1]=$2
}
ARGIND==2 {
    coverage=$2;
    split($NF, reads, ","); 
    for (i=3; i<=2+coverage; i++) {
        rname = reads[i-2];
        if ($1==POS) focal_reads[rname]=1;
        base[rname][$1]=$i;
    } 
}
END {
    for(rname in focal_reads) {
        if(length(base[rname]) == 1) continue;
        focal_phase = (base[rname][POS] == "." || base[rname][POS] == ",") ? 0x1 : 0x2;
        for(pos in base[rname]) {
            if (pos == POS) continue;
            phaseable_positions[pos]=1;
            current_phase = (base[rname][pos] == "." || base[rname][pos] == ",") ? 0x1 : 0x2;
            relative_phase = or(focal_phase, current_phase)
            if (and(relative_phase, 0x3) == 0x3) 
                opposite_phase[pos]+= 1;
            else
                same_phase[pos]+=1;
        }
    }
    for (variant in variant2pos) {
        pos = variant2pos[variant];
        if (pos == POS ) {
            GT="0|1";
            PS=POS;
        }
        else if(same_phase[pos]+0 > 0 && opposite_phase[pos]+0 > 0 ) {
            GT="0/1";
            PS=".";
        } else if (same_phase[pos]+0 > 0) {
            GT="0|1";
            PS=POS;
        } else if (opposite_phase[pos]+0 > 0) {
            GT="1|0";
            PS=POS;
        } else {
            GT="0/1";
            PS=".";
        }
        print ID, variant, GT, PS, same_phase[pos]+0, opposite_phase[pos]+0;
    }
}
' \
    <(bcftools query -f '[%CHROM:%POS:%REF:%ALT\t%POS\n]' $TMP_DIR/$PREFIX.input.vcf ) \
    -  \
    | sort -k2,2V

rm -f $TMP_DIR/$PREFIX.input.vcf

# set on input:
# CHR

set -euo pipefail
set +e 
sudo apt-get update 2> /dev/null 1>&2
set -e
sudo apt --yes install parallel 2> /dev/null 1>&2

export TMP_DIR=$HOME/tmp
mkdir -p $TMP_DIR


localize () {
    FILENAME="$TMP_DIR/`basename $1`"
    [ ! -f $FILENAME ] && cp $1 $TMP_DIR
    echo $FILENAME
}

export REFERENCE=`localize /mnt/project/lohdata/resources/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa`
export REFERENCE_INDEX=`localize /mnt/project/lohdata/resources/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai`
export RARE_VAR_VCF=`localize /mnt/project/lohdata/resources/burden_masks/vcfs/$CHR.LoF.scaffold.phased.AS.bcf`
export RARE_VAR_VCF_CSI=`localize /mnt/project/lohdata/resources/burden_masks/vcfs/$CHR.LoF.scaffold.phased.AS.bcf.csi`
export BURDEN_MASKS_MASKS=`localize /mnt/project/lohdata/resources/burden_masks/mask_definitions/$CHR.gnomAD.LoF.missense.CNV.masks`
export BURDEN_MASKS_ANNOTATIONS=`localize /mnt/project/lohdata/resources/burden_masks/mask_definitions/$CHR.gnomAD.LoF.missense.CNV.annotations.txt.gz`
export PROTEIN_CODING_VARIANTS=`localize /mnt/project/lohdata/resources/burden_masks/mask_definitions/gnomAD.protein_coding.high_moderate_vep.$CHR.txt.gz`
export BURDEN_MASKS_DIR=$TMP_DIR
export ALLELIC_SHIFT_VCF=`localize /mnt/project/lohdata/david/mCAs_WGS/calls/allelic_shift/WGS_500k.$CHR.allelic_shift.bcf`
export ALLELIC_SHIFT_VCF_CSI=`localize /mnt/project/lohdata/david/mCAs_WGS/calls/allelic_shift/WGS_500k.$CHR.allelic_shift.bcf.csi`
export CALLS_FILE=`localize /mnt/project/lohdata/david/mCAs_WGS/calls/WGS_500k.calls.txt`

OUT_PFX=$OUT_DIR/ID.variant.mask.read_phase.statistical_phase.$CHR
RESULTS_DIR=/mnt/project/lohdata/david/mCAs_WGS/cisGWAS/results
SUFFIX=LoF.protein_coding.CN-LOH.GWAS.txt.gz

zcat $RESULTS_DIR/$CHR.$SUFFIX \
    | awk -v OFS='\t' -v CHR=$CHR '
$6>0 && $1==CHR && $5=="LoF.all_transcripts.missense8.0.001" {
    split($3,a,"."); 
    gene=a[1]; 
    if(gene in seen) next; 
    seen[gene]=1; 
    mask[gene]=a[2]; 
    for(i in a) {
        if(i<3) continue; 
        if(a[i]=="singleton") {
            AF[gene]="singleton";
            break;
        }
        if(a[i]==0) {
            AF[gene]=a[i]"."a[i+1];
            break
        }; 
        mask[gene]=mask[gene]"."a[i]
    }; 
    if(!(gene in AF)) next;
    counter+=1;
    print $1,$2,gene,mask[gene],AF[gene];
} ' \
    | sort -k2,2n \
    | parallel --colsep='\t' --retries 4 --halt-on-error 2 "CHR={1} GENE={3} MASK={4} AF={5} bash refine_phase.sh" \
    | tee $OUT_PFX.txt  \
    | awk -v OFS='\t' '
BEGIN {
    print "mask","carrier_overrep","carrier_underrep","missing_phase","inconsistent_phase";
}
{
    masks[$3]=1; 
    phase=and(or($4,$5), 0x3); 
    if (phase==0x1) overrep[$3]++; 
    else if (phase==0x2) underrep[$3]++; 
    else if (phase==0) missing_phase[$3]++; 
    else inconsistent_phase[$3]++;
} 
END { 
    for (mask in masks) print mask,overrep[mask]+0,underrep[mask]+0,missing_phase[mask]+0,inconsistent_phase[mask]+0
}' \
    | gzip > $OUT_DIR/$CHR.burden.allelic_shift.GWAS.txt.gz

gzip $OUT_PFX.txt

# require gnu-parallel to be installed
# set on input
# CHR
# GENE
# MASK
# AF
# REFERENCE

set -euo pipefail

TMP_DIR=$HOME/tmp
mkdir -p $TMP_DIR

PHASING_SCRIPT=${PHASING_SCRIPT:-read_based_phasing.sh}
echo "`date`: performing read based phasing with $PHASING_SCRIPT..." 1>&2 
[[ $PHASING_SCRIPT == "run_whatshap_1.sh" ]] && pip3 install whatshap 1>&2

GENE_MASK_AF=$GENE.$MASK.$AF
RARE_VAR_VCF=${RARE_VAR_VCF:-/mnt/project/lohdata/resources/burden_masks/vcfs/$CHR.LoF.scaffold.phased.AS.bcf}
BURDEN_MASKS_DIR=${BURDEN_MASKS_DIR:-/mnt/project/lohdata/resources/burden_masks/mask_definitions}
ALLELIC_SHIFT_VCF=${ALLELIC_SHIFT_VCF:-/mnt/project/lohdata/david/mCAs_WGS/calls/allelic_shift/WGS_500k.$CHR.allelic_shift.bcf}
CALLS_FILE=${CALLS_FILE:-/mnt/project/lohdata/david/mCAs_WGS/calls/WGS_500k.calls.txt}
PREFIX=$GENE_MASK_AF
EVENT_TYPE=${EVENT_TYPE:-CN-LOH}

# extract a region that contains all the protein coding variants for the gene and nearby common variants
REGION=`zcat $BURDEN_MASKS_DIR/gnomAD.protein_coding.high_moderate_vep.$CHR.txt.gz \
    | awk -v GENE=$GENE '
BEGIN {
    minpos=1e9;
    maxpos=0;
}
$10==GENE {
    chrom=$1;
    minpos=($2<minpos) ? $2: minpos; 
    maxpos=($2>maxpos) ? $2 : maxpos
} 
END {
    print chrom":"minpos-10000"-"maxpos+10000
}'`


# keep only individuals with a CN-LOH overlapping the region of interest
awk -v REGION=$REGION -v CHR=$CHR -v EVENT_TYPE=$EVENT_TYPE '
FNR==1 {
    for(i=1; i<=NF; i++) f[$i]=i;
    split(REGION, a, ":");
    split(a[2], b, "-");
    bpStart=b[1];
    bpEnd=b[2];
}
ARGIND==1 {keep[$1]=1}
ARGIND==2 && $f["chr"] == CHR && $f["type"]==EVENT_TYPE {   
    if($f["p"] == "T") $f["bpStart"] = 1;
    if($f["q"] == "T") $f["bpEnd"] = 1e9;
    if($f["bpStart"] < bpEnd && $f["bpEnd"] > bpStart) 
        if($f["ID"] in keep) 
            print $f["ID"];
}' \
    <(bcftools query -l $RARE_VAR_VCF) \
    $CALLS_FILE \
    | sort -u \
    > $TMP_DIR/$PREFIX.sample_list.txt 

if [[ `wc -l $TMP_DIR/$PREFIX.sample_list.txt` == 0 ]]; then
    echo "error: no individuals with mCAs overlapping $GENE" 2>&1 
    exit 1
fi

# extract LoF variant for all individuals
bcftools query -S $TMP_DIR/$PREFIX.sample_list.txt -r $REGION -f '[%SAMPLE\t%ID\t%GT\n]' $RARE_VAR_VCF \
    | grep -w $GENE_MASK_AF \
    | awk '($3=="1|0" || $3=="0|1" || $3=="0/1" || $3=="1/0") {print $1}' \
    | parallel "REGION=$REGION GENE=$GENE MASK=$MASK BURDEN_MASKS_DIR=$BURDEN_MASKS_DIR ID={} bash extract_LoF_variant.sh" \
    | tee $TMP_DIR/$PREFIX.LoF_variant.txt \
    | sort -k1,1 -u \
    | parallel --joblog=/dev/stderr --colsep='\t' "ID={1} CHR={2} POS={3} REFERENCE=$REFERENCE bash $PHASING_SCRIPT" \
    > $TMP_DIR/$PREFIX.phase.txt

# align LoF phase with mCA phase
# output code: 
# 0  -- missing variant (CNV perhaps)
# 1  -- LoF made hom
# 2  -- LoF removed
# 4  -- missing phase
# 8  -- unable to merge phase into scaffold
# 16 -- inconsistent phase
awk -v OFS='\t' -v GENE_MASK_AF=$GENE_MASK_AF '
function which_hap_overrep(AS, GT) {
    if((AS==-1 && GT=="0|1") || (AS==1 && GT=="1|0")) {
        overrep_hap=1;
    } else if((AS==-1 && GT=="1|0") || (AS==1 && GT=="0|1")) {
        overrep_hap=0;
    } else {
        overrep_hap=-9;
    }
    return overrep_hap;
}
function get_phase_stats(hap_phase, GT) {
    if ((hap_phase==1 && GT=="0|1") || (hap_phase==0 && GT=="1|0")) res = 0x1;
    else if ((hap_phase==1 && GT=="1|0") || (hap_phase==0 && GT=="0|1")) res = 0x2;
    else res = 0x10;
    return res;
}
ARGIND==1 && ($3=="1|0" || $3=="0|1") && ($4==1 || $4==-1) {
    GT_ref[$1,$2]=$3;
    AS_ref[$1,$2]=$4;
    overrep_hap=which_hap_overrep($4,$3);
    if(($1 in hap_phase) && hap_phase[$1]!=overrep_hap) hap_phase[$1]=-9;
    else hap_phase[$1]=overrep_hap;
} 
ARGIND==2 {   
    ID_to_variant[$1]=$2":"$3":"$4":"$5
}
ARGIND==3 && $4!="." {
    phase_group_size[$1,$4]+=1;
    if (!($1,$2) in AS_ref) {
        next;
    }
    GT=$3;
    AS=AS_ref[$1,$2];
    overrep_hap = which_hap_overrep(AS, GT);

    if((($1,$4) in phase_group) && phase_group[$1,$4] != overrep_hap) {
        phase_group[$1,$4]=-9;
    } else {
        phase_group[$1,$4]=overrep_hap;
    }
}
ARGIND==4 {
    if(ID_to_variant[$1] != $2) next;
    if(!(($1,$4) in phase_group)) read_phasing[$1] = (phase_group_size[$1,$4]+0>1) ? 0x8 : 0x4;
    else if (phase_group[$1,$4] == -9) read_phasing[$1] = 0x10;
    else read_phasing[$1] = get_phase_stats(phase_group[$1,$4], $3);
}
ARGIND==5 {
    if ($2!=GENE_MASK_AF) next;
    if(!($1 in hap_phase) || ($4<0.8)) {
        statistical_phasing[$1] = 0x4;
    }
    else if (hap_phase[$1] == -9) statistical_phasing[$1] = 0x10;
    else statistical_phasing[$1] = get_phase_stats(hap_phase[$1], $3);
    ID_to_variant[$1]=($1 in ID_to_variant) ? ID_to_variant[$1] : ".";
}
END {
    for(ID in ID_to_variant) {
        print ID,ID_to_variant[ID],GENE_MASK_AF,read_phasing[ID]+0,statistical_phasing[ID]+0;
    }
}
' \
    <( bcftools view -S $TMP_DIR/$PREFIX.sample_list.txt $ALLELIC_SHIFT_VCF -r $REGION | bcftools query -f '[%SAMPLE\t%CHROM:%POS:%REF:%ALT\t%GT\t%AS\n]' ) \
    $TMP_DIR/$PREFIX.LoF_variant.txt \
    $TMP_DIR/$PREFIX.phase.txt \
    $TMP_DIR/$PREFIX.phase.txt \
    <( bcftools view -S $TMP_DIR/$PREFIX.sample_list.txt $RARE_VAR_VCF -r $REGION | egrep "#|$GENE_MASK_AF"  | bcftools query -f '[%SAMPLE\t%ID\t%GT\t%PP\n]' | awk '$3=="0|1" || $3=="1|0"' ) 

rm -f $TMP_DIR/$PREFIX.{sample_list,phase,LoF_variant}.txt 

# set on input
# BFILE

set -euo pipefail

TMP_DIR=$HOME/tmp
mkdir -p $TMP_DIR

INPUT_DIR=/mnt/project/lohdata/resources/burden_masks/mask_definitions

plink2 \
    --bfile $BFILE \
    --freq \
    --remove /mnt/project/lohdata/resources/BOLT-LMM/remove.nonEUR.FID_IID.40709.txt \
    --out $TMP_DIR/$CHR.ukb_AF

zcat $INPUT_DIR/gnomAD.protein_coding.high_moderate_vep.$CHR.txt.gz \
    | sed 's/^chr//' \
    | awk -v OUT_DIR=$OUT_DIR -v CHR=$CHR '
ARGIND==1 {
    variant=$2;
    ukb_AF[variant]=$5;
}
ARGIND==2 && ($14 == "LoF" || index($11, "missense_variant") != 0) {
    total_variants+=1;
    variant = $1":"$2":"$4":"$5;
    gene = $10;
    consequence = $11;
    primateAI = $15;
    non_ukb_nfe_AF = ($8!=0) ? ($9+1)/$8 : 0;
    if (sets[gene]=="") {
        chr_pos[gene] = $1"\t"$2;
        sets[gene] = variant;
    }
    else {
        sets[gene] = sets[gene]","variant;
    }
    if (index(consequence, "missense_variant") != 0) {
        effect="missense"int(primateAI*10);
    } 
    else {
        effect="LoF";
    }
    if ($13!="MANE_SELECT") {
        effect=effect".NON_MANE_SELECT"
    }
    if ((non_ukb_nfe_AF*10 < ukb_AF[variant]+0) && (ukb_AF[variant]+0 > 1e-5)) {
        printf "skipping %s (UKB_AF=%s, gnomAD_AF=%s)\n", variant, ukb_AF[variant], non_ukb_nfe_AF > "/dev/stderr"
        skipped_variants+=1;
        next;
    }
    print variant,gene,effect;
    cnv[gene] = "chr"$1"."gene".pLoF_highConf";
}
END {
    printf "%s/%s variants skipped\n", skipped_variants,total_variants > "/dev/stderr"
    for (gene in sets) {
        print cnv[gene],gene,"LoF_CNV";
        print gene"\t"chr_pos[gene]"\t"sets[gene]","cnv[gene] | "sort -k2,2n -k3,3n | gzip > " OUT_DIR"/"CHR".gnomAD.LoF.missense.CNV.sets.txt.gz"
    }
}' $TMP_DIR/$CHR.ukb_AF.afreq - \
    | gzip > $OUT_DIR/$CHR.gnomAD.LoF.missense.CNV.annotations.txt.gz

awk -v CHR=$CHR '
function push(arr, x) {
    arr[length(arr)+1] = x
}
function print_mask(MANE_SELECT, CNV, MISSENSE_SCORE) {
    mask_name="LoF"
    if (MANE_SELECT==0) mask_name=mask_name".all_transcripts"
    if (CNV==1) mask_name=mask_name".CNV"
    if (MISSENSE_SCORE<=8) mask_name=mask_name".missense"MISSENSE_SCORE
    printf "%s\t", mask_name;


    arr[1] = "LoF"
    if (MISSENSE_SCORE <= 8) {
        for(i=9; i>= MISSENSE_SCORE; i--) {
            push(arr, "missense"i);
        }
    }
    
    if(MANE_SELECT==0){
        num_masks=length(arr);
        for(i=1; i<=num_masks; i++) {
            push(arr, arr[i]".NON_MANE_SELECT");
        }
    }

    if (CNV==1) {
        push(arr, "LoF_CNV");
    }
    for (i=1; i<=length(arr); i++) {
        if(i!=1) printf ","
        printf "%s", arr[i]
    }
    printf "\n"
}
BEGIN {
    for (MISSENSE_SCORE=9; MISSENSE_SCORE>=6; MISSENSE_SCORE--) {
        for (MANE_SELECT=0; MANE_SELECT<=1; MANE_SELECT++) {
            for (CNV=0; CNV<=1; CNV++) {
                if (CHR=="chrX" && CNV==1) continue;
                print_mask(MANE_SELECT, CNV, MISSENSE_SCORE);
                delete arr;
            }
        }
    }
}' > $OUT_DIR/$CHR.gnomAD.LoF.missense.CNV.masks

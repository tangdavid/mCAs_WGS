# set on input:
# CHR

set -euo pipefail

TMP_DIR=$HOME/tmp
mkdir -p $TMP_DIR

echo "`date`: subsetting to individuals with calls..." 1>&2
CALL_FILE=/mnt/project/lohdata/david/mCAs_WGS/calls/WGS_500k.calls.txt
awk -v chrom=$CHR '
BEGIN {call_file_num = 5}
ARGIND<call_file_num {
    exclude[$1]=1
}
ARGIND==call_file_num && FNR==1 {
    for(i=1; i<= NF; i++) f[$i]=i
} 
ARGIND==call_file_num && FNR>1 && !($1 in exclude) {
    if( $(f["chr"])==chrom && $(f["type"])=="CN-LOH" && $(f["length"])>1e6) {
        print $1
    }
}' \
    /mnt/project/lohdata/resources/BOLT-LMM/remove.nonEUR.FID_IID.40709.txt \
    /mnt/project/lohdata/resources/BOLT-LMM/w40709_CURRENT.FID_IID.txt \
    <( CHR=$CHR bash ../cisGWAS/remove_rare_variant_carriers.sh ) \
    <( cut -f1 /mnt/project/lohdata/david/mCAs_WGS/calls/WGS_500k.calls.txt | sort | uniq -c | awk '$1>1 {print $2}' ) \
    $CALL_FILE \
    | bcftools view /mnt/project/lohdata/david/mCAs_WGS/polygenic_drive/snp_array_bcf/$CHR.snp_array.AS.bcf \
        -S - \
        -Ou -o $TMP_DIR/$CHR.snp_array.single_event.EUR.bcf 

function run_single_trait() {
    TRAIT=$1
    BCF_FILE=$2
    ../../bin/mCAs_WGS compute-polygenic-drive \
        --prs /mnt/project/lohdata/david/mCAs_WGS/polygenic_drive/predbetas/bolt_460K_selfRepWhite.$TRAIT.predbetas.txt.gz \
        < $BCF_FILE
}

export -f run_single_trait

ls /mnt/project/lohdata/david/mCAs_WGS/polygenic_drive/predbetas/ \
    | cut -f2 -d'.' \
    | parallel --joblog=/dev/stderr "run_single_trait {} $TMP_DIR/$CHR.snp_array.single_event.EUR.bcf > $TMP_DIR/$CHR.{}.differential_PRS.txt" 


awk -v OFS='\t' '
BEGIN {
    print "ID","event_arm","trait","differential_PRS"
}
ARGIND==1 && FNR==1 {for(i=1; i<= NF; i++) f[$i]=i} F
ARGIND==1 && FNR>1 {
    event_arm=$(f["chr"]); 
    if($(f["p"])!="N") event_arm=event_arm"p"; 
    if($(f["q"])!="N") event_arm=event_arm"q"; 
    label[$(f["ID"])]=event_arm
}
ARGIND>1 && FNR > 1 {
    split(FILENAME, a, "."); 
    trait=a[2]; 
    traits[trait]=1; 
    scores[$1][trait] = $2
} 
END {
    for(trait in traits) 
        for (ID in scores) 
            print ID, label[ID], trait, scores[ID][trait]
}' \
    $CALL_FILE \
    $TMP_DIR/$CHR.*.differential_PRS.txt \
    > $OUT_DIR/$CHR.differential_PRS.txt 

rm -f $TMP_DIR/$CHR.*.differential_PRS.txt
rm -f $TMP_DIR/$CHR.snp_array.single_event.EUR.bcf


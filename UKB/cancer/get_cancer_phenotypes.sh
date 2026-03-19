
set -euo pipefail

TMP_DIR=$HOME/tmp
mkdir -p $TMP_DIR

dx extract_dataset /app40709_20250401062123.dataset --fields participant.eid,participant.p3166_i0_a0 -o - \
    | awk -F, '$2!="" {split($2,a," "); print $1","a[1]}' > $TMP_DIR/sample_date.txt 
dx extract_dataset /app40709_20250401062123.dataset --fields `echo participant.eid participant.p40005_i{0..21} | tr ' ' ','` -o - \
    > $TMP_DIR/cancer_date.txt 
dx extract_dataset /app40709_20250401062123.dataset --fields `echo participant.eid participant.p40012_i{0..21} | tr ' ' ','` -o - \
    > $TMP_DIR/behavior.txt 
dx extract_dataset /app40709_20250401062123.dataset --fields `echo participant.eid participant.p40011_i{0..21} | tr ' ' ','` -o - \
    > $TMP_DIR/histology.txt 

awk -F, '
BEGIN {
    OFS="\t"
    print "ID","histology","behavior","collection_date","cancer_date"
}
ARGIND==1 {keep[$1]++} 
ARGIND==2 {blood_date[$1]=$2}
ARGIND==3 {for (i=2;i<=NF;i++) if($i!="") cancer_date[$1][i]=$i}
ARGIND==4 {for(i=2;i<=NF;i++) if($i!="") behavior[$1][i]=$i}
ARGIND==5 && FNR>1 {
    for (i=2;i<=NF;i++) {
        if ($i=="") continue
        print $1,$i,behavior[$1][i],($1 in blood_date) ? blood_date[$1] : "NA",cancer_date[$1][i]
    } 
}'  \
    /mnt/project/lohdata/david/mCAs_WGS/sample_data/WGS_500k.included_IDs.txt \
    $TMP_DIR/sample_date.txt \
    $TMP_DIR/cancer_date.txt \
    $TMP_DIR/behavior.txt \
    $TMP_DIR/histology.txt \
    > ID.histology.behavior.sample_date.cancer_date.txt

set -euo pipefail

TMP_DIR=$HOME/tmp
mkdir -p $TMP_DIR
BREAKPOINT_FILE=$1

cat $BREAKPOINT_FILE \
    | awk -v OFS='\t' '
$1!="ID" {
    name=$1
    print $3,$4-50,$4,name;
    print $3,$4,$4+50,name;
    print $7,$8-50,$8,name;
    print $7,$8,$8+50,name;
}' \
    | tr ':-' '_' \
    | bedtools getfasta -fi ~/tmp/GRCh38_full_analysis_set_plus_decoy_hla.fa -bed - -name \
    | awk 'index($1,">") {$1=substr($1,1,8)} {print}'

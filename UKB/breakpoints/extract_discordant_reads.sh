set -euo pipefail
set +e 
sudo apt-get update 1>&2
set -e
sudo apt --yes install parallel 1>&2

EVENTS_FILE=$1
TMP_DIR=$HOME/tmp
mkdir -p $TMP_DIR
cp /mnt/project/lohdata/resources/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa{,.fai} $TMP_DIR

awk '
function max(x, y) {
    return (x > y) ? x:y
}
$1!="ID" {
    arr[$1]=arr[$1]","$2":"max($3-50e3,1)"-"$3+50e3","$2":"max($4-50e3,1)"-"$4+50e3; 
    seen[$1]+=1
} 
END {
    for (ID in arr) {
        if(seen[ID]<20) print ID,substr(arr[ID],2) 
    }
}' $EVENTS_FILE \
    | parallel -j16 --colsep=' ' --halt-on-error 2 --joblog=/dev/stderr REF_DIR=$TMP_DIR bash reconstruct_complex_sv.sh {1} {2} 

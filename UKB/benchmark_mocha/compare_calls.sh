set -euo pipefail

BENCHMARK=$1

TMP_DIR=$HOME/tmp
mkdir -p $TMP_DIR

rm -f mocha_calls.txt deduped_ads.txt

awk 'ARGIND==1{remove[$1]++} ARGIND==2 && !($1 in remove) && $3<=2 && $2!~"WGS7|WGS8|WGS9" {print}' \
    /mnt/project/lohdata/resources/BOLT-LMM/w40709_CURRENT.FID_IID.txt \
    /mnt/project/lohdata/resources/WGS_depth/sample_batches/ID_40709.WGS_group.batch.txt \
    > $TMP_DIR/mocha.benchmark.ind.txt

awk 'ARGIND==1 {keep[$1]++} ARGIND==2 && ($1 in keep) {print $1,$3}' \
    $TMP_DIR/mocha.benchmark.ind.txt \
    /mnt/project/lohdata/david/mCAs_WGS/calls/WGS_500k.calls.txt  \
    | sort  \
    > deduped_ads.txt 

awk 'ARGIND==1 {keep[$1]++} ARGIND==2 && ($1 in keep) && $3!="chrom" && $3!="chrX" {print $1,$3}' \
    $TMP_DIR/mocha.benchmark.ind.txt \
    <(zcat /mnt/project/lohdata/david/mCAs_WGS/benchmark_mocha/WGS_10k.mocha.$BENCHMARK.calls.txt.gz) \
    | sort  \
    > mocha_calls.txt 

set +e
diff mocha_calls.txt deduped_ads.txt \
    | awk -v BENCHMARK=$BENCHMARK -v OFS='\t' '
ARGIND==1 && $1==">" {new++} 
ARGIND==1 && $1=="<" {old++} 
ARGIND==2 {mocha_calls++}
ARGIND==3 {deduped_ads++}
END {
    print "calls from MoChA with "BENCHMARK" settings:", mocha_calls
    print "calls from full WGS pipeline:", deduped_ads
    print "exclusive to full WGS pipeline:",new
    print "exclusive to MoChA pipeline:",old
}' \
    - \
    mocha_calls.txt \
    deduped_ads.txt 
set -e
rm mocha_calls.txt deduped_ads.txt

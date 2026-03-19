set -euo pipefail

BUFFER=${1:-2}

DISC_BREAKPOINT_FILE=/mnt/project/lohdata/david/mCAs_WGS/breakpoints/WGS_500k.discordant_breakpoints.txt.gz

zcat $DISC_BREAKPOINT_FILE \
    | awk '$1!="ID" && ($7=="INNER" && $11=="LOSS" || $7=="OUTER" && $11=="GAIN") {print $1,$2":"$5"-"$6}' \
    | parallel -j16 --colsep=' ' --joblog=/dev/stderr bash get_split_reads.sh {1} {2} $BUFFER \
    > $OUT_DIR/WGS_500k.calls.all_split_reads.txt

bash compute_microhomology.sh $OUT_DIR/WGS_500k.calls.all_split_reads.txt \
    | awk -v OFS='\t' -F'\t' '
function abs(x) {return x>0 ? x : -x} 
function max(x, y) {return x>y ? x : y} 
BEGIN {
    print "ID","chr","bpStart","bpEnd","orientation","leftBreakpoint","rightBreakpoint","microhomology","microhomologySeq","numberReads","potentialBreakpoints";
}
ARGIND==1 {
    callMap[$1"_"$2":"$5"-"$6] = $1 OFS $2 OFS $3 OFS $4 OFS $7
}
ARGIND==2 {
    ID=$1;

    call = callMap[ID"_"$2];

    leftBreakpoint=$4;
    rightBreakpoint=$8;
    microhomology=$12;
    microhomologySeq=$13;

    exact_breakpoint=call OFS leftBreakpoint OFS rightBreakpoint OFS microhomology OFS microhomologySeq
    seen[call]+=(exact_breakpoint in num_reads) ? 0 : 1;
    num_reads[exact_breakpoint]+=1; 
    exact_breakpoint_to_call[exact_breakpoint]=call;
}
END {
    for (exact_breakpoint in num_reads) {
        print(exact_breakpoint, num_reads[exact_breakpoint], seen[exact_breakpoint_to_call[exact_breakpoint]]) | "sort" 
    }
}
' <(zcat $DISC_BREAKPOINT_FILE) - \
    > $OUT_DIR/WGS_500k.split_breakpoints.txt

gzip -f $OUT_DIR/WGS_500k.split_breakpoints.txt
gzip -f $OUT_DIR/WGS_500k.calls.all_split_reads.txt

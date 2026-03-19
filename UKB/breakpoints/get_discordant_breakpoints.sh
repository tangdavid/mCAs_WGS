set -euo pipefail

DISCORDANT_READS_FILE=$1

cat $DISCORDANT_READS_FILE \
    | awk -v OFS='\t' '
$2==$4 {
    if($3 < $5) {
        start=$3; 
        end=$5;
        left_forward=$6;
        right_forward=$7;
    } 
    if($5<$3) {
        start=$5; 
        end=$3;
        left_forward=$7;
        right_forward=$6;
    } 
    print $1,$2,start,end,left_forward,right_forward,$8
}' \
    | bash merge_adjacent_bins.sh \
    | awk '
BEGIN {OFS="\t"; print "ID","chr","bpStart","bpEnd","breakpoints","leftForward","rightForward","totalReads"} 
(($5==$7 && $6==0) || ($6==$7 && $5==0) || ($5==0 && $6==0)  || ($5==$6 && $6==$7)){
    print $1,$2,$3,$4,$2":"$3"-"$4,$5,$6,$7
}' \
    | bash merge_breakpoints_with_calls.sh  /mnt/project/lohdata/david/mCAs_WGS/calls/WGS_500k.calls.txt - \
    | tr ':-' '\t' \
    | awk -v OFS='\t' '
function abs(x) {return x>0 ? x : -x} 
function check_orientation(orientation, cn) {
    return (orientation=="INNER" && type[ID_region]=="LOSS") || (orientation=="OUTER" && type[ID_region]=="GAIN")
}
BEGIN {
    print "ID","chr","bpStart","bpEnd","leftBreakpoint","rightBreakpoint","orientation","totalReads","combinedError","potentialBreakpoints","type"
}
ARGIND==1 && FNR==1 { for(i=1;i<=NF;i++) f[$i]=i}
ARGIND==1 && FNR>1 {
    ID_region = $f["ID"] OFS $f["chr"] OFS $f["bpStart"] OFS $f["bpEnd"]
    type[ID_region] = $f["type"]
}
ARGIND==2 && FNR>1 {
    region=$2 OFS $3 OFS $4; 
    current_dist=abs($3-$6)+abs($4-$7); 
    ID_region=$1 OFS region;
    seen[ID_region] ++;
    if ($8==$10 && $9==0) {
        orientation="INNER"
    }
    else if ($9==$10 && $8==0) {
        orientation="OUTER"
    }
    else {
        orientation="UNKNOWN"
    }
    proper_pairing[ID_region] += (check_orientation(orientation, type[ID_region])) ? 1 : 0;
    if(current_dist < best_dist[ID_region] || seen[ID_region]==1 || proper_pairing[ID_region]==1) {
        if(seen[ID_region]>1 && !check_orientation(orientation, type[ID_region])) next;
        breakpoint[ID_region]=$6 OFS $7 OFS orientation OFS $10; 
        best_dist[ID_region]=current_dist
        best_orientation=orientation
    }
} 
END {
    for (ID_region in breakpoint) 
        print ID_region,breakpoint[ID_region],best_dist[ID_region],proper_pairing[ID_region],type[ID_region];
}' /mnt/project/lohdata/david/mCAs_WGS/calls/WGS_500k.calls.txt -

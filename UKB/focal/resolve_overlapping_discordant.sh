# set on input:
# Z_SCORE
set -euo pipefail 

Z_SCORE=${Z_SCORE:-3}

awk -v OFS='\t' '
BEGIN {
    print "ID","chr","bpStart","bpEnd","OBSinclCNV","EXPinclCNV","orientation","region"
}
ARGIND==1 {
    keep[$1]++;
}
ARGIND==2 && FNR==1 {
    for(i=1; i<=NF; i++) f[$i]=i
} 
ARGIND==2 && FNR>1 && $1 in keep {
    region=$(f["region"]);
    split(region,a,":");
    split(a[2],b,"-");
    windowStart=b[1];
    windowEnd=b[2];
    print $(f["ID"]),$(f["chr"]),$(f["bpStart"]),$(f["bpEnd"]),$(f["OBSreads"]),$(f["EXPreads"]),"INNER",$(f["region"]);
    print $(f["ID"]),$(f["chr"]),windowStart-5e6,windowEnd+5e6,$(f["OBSsink"]),$(f["EXPsink"]),"SINK",$(f["region"]);
}' <(zcat /mnt/project/lohdata/david/mCAs_WGS/del13q/WGS_500k.del13q.calls.depth.txt.gz) - \
    | sort -k1,1 -k2,2 -k3,3 -k4,4 -un \
    | python3 resolve_overlapping_discordant.py  --zscore $Z_SCORE

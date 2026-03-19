# set on input:
# BIN_SIZE

set -euo pipefail 

awk -v OFS='\t' -v BIN_SIZE=$BIN_SIZE '
BEGIN{
    BIN_SIZE=BIN_SIZE*1e3
}
NR==1 {
    for(i=1; i<=NF; i++) {
        f[$i]=i;
    }
    print $0
} 
NR>1 { 
    ID=$(f["ID"]);
    region=$(f["region"]);
    start=int($(f["bpStart"])/BIN_SIZE)*BIN_SIZE; 

    if(start != prevStart && NR>2) {
        if(ID==prevID && region==prevRegion)
            print ID,region,prevStart,start,obs,expected,expected_noPC,(expected!=0)?(obs/expected):"-nan",0; 
        obs=0; 
        expected=0;
        expected_noPC=0;
    } 
    if ($(f["inCNV"])==0) {
        obs+=$(f["OBSreads"]); 
        expected+=$(f["EXPreads"]); 
        expected_noPC+=$(f["EXPreadsNoPCadj"]);
    }
    prevStart=start;
    prevID=ID;
    prevRegion=region;

}
END {
    print ID,region,prevStart,prevStart+BIN_SIZE,obs,expected,expected_noPC,(expected!=0)?(obs/expected):-nan,0; 
}
' /dev/stdin

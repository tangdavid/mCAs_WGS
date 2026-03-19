set -euo pipefail

DISCORDANT_READ_FILE=$1
PROFILE_FILE=$2

awk -v OFS='\t' '
BEGIN{ buffer=100e3}
ARGIND==1 {
    split($2,a,":"); 
    split(a[2],b,"-"); 
    discordant[$1][$2]=1; 
    chrom[$1,$2]=a[1]; 
    bpStart[$1,$2]=b[1]; 
    bpEnd[$1,$2]=b[2];

    meta[$1,$2]=$3 OFS $5;
} 

ARGIND==2 && FNR==1 { for (i=1; i<=NF; i++) f[$i]=i; }
ARGIND==2 && ($1 in discordant) {
    ID=$(f["ID"]); 
    split($(f["region"]),a,":"); 
    currentChrom=a[1]; 
    for (region in discordant[ID]) {
        if(chrom[ID,region] != currentChrom) continue;
        if($(f["bpStart"])>=bpStart[ID,region] && $(f["bpEnd"])<=bpEnd[ID,region]) {
            obs_region[ID,region]+=$(f["OBSreads"]);
            exp_region[ID,region]+=$(f["EXPreads"]);
        }
        if($(f["bpStart"])>=bpStart[ID,region]-buffer && $(f["bpEnd"])<=bpEnd[ID,region]+buffer) {
            obs_outer[ID,region]+=$(f["OBSreads"]);
            exp_outer[ID,region]+=$(f["EXPreads"]);
        }
        if($(f["bpStart"])>=bpStart[ID,region]+buffer && $(f["bpEnd"])<=bpEnd[ID,region]-buffer) {
            obs_inner[ID,region]+=$(f["OBSreads"]);
            exp_inner[ID,region]+=$(f["EXPreads"]);
        }
        if($(f["bpStart"])>=bpStart[ID,region]-5e6 && $(f["bpEnd"])<=bpEnd[ID,region]+5e6) {
            obs_sink[ID,region]+=$(f["OBSreads"]);
            exp_sink[ID,region]+=$(f["EXPreads"]);
        }
    }
}
END {
    print "ID","chr","bpStart","bpEnd","OBSreads","EXPreads","OBSflank","EXPflank","OBSinner","EXPinner","OBSsink","EXPsink","totalReads","region"
    for (ID in discordant) {
        for (region in discordant[ID]) {
            print \
                ID,\
                chrom[ID,region],\
                bpStart[ID,region],\
                bpEnd[ID,region],\
                obs_region[ID,region],\
                exp_region[ID,region],\
                obs_outer[ID,region]-obs_region[ID,region],\
                exp_outer[ID,region]-exp_region[ID,region],\
                obs_region[ID,region]-obs_inner[ID,region],\
                exp_region[ID,region]-exp_inner[ID,region],\
                obs_sink[ID,region],\
                exp_sink[ID,region],\
                meta[ID,region]
        }
    }
}
' \
    $DISCORDANT_READ_FILE $PROFILE_FILE

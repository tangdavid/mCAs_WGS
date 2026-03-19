# set on input:
# CALLS_FILE
# DEPTH_PROFILES

set -euo pipefail


awk -F'\t' -v OFS='\t' '
ARGIND==1 {cen[$1]=$3} 
ARGIND==2 && FNR==1 {for(i=1; i<=NF; i++) f2[$i]=i}
ARGIND==4 && FNR==1 {for(i=1; i<=NF; i++) f4[$i]=i; print $0}
ARGIND==2 {
    if ($f2["p"]=="T" && $f2["q"]=="T") {
        observed[$1][$3]["p"]=0; 
        observed[$1][$3]["q"]=0; 

        expected[$1][$3]["p"]=0;
        expected[$1][$3]["q"]=0;
    }
}
ARGIND==3 && $1 in observed && $2 in observed[$1] {
    arm=($4<cen[$2]) ? "p" : "q"
    observed[$1][$2][arm] += $5;
    expected[$1][$2][arm] += $6;
}
ARGIND==4 && FNR>1 {
    skip=0;
    ID=$f4["ID"]
    chr=$f4["chr"]
    if(ID in observed && chr in observed[ID]) {
        obs_p = observed[ID][chr]["p"];
        obs_q = observed[ID][chr]["q"];
        exp_p = expected[ID][chr]["p"];
        exp_q = expected[ID][chr]["q"];

        z_p = (obs_p - exp_p)/sqrt(obs_p);
        z_q = (obs_q - exp_q)/sqrt(obs_q);

        if (z_p * z_q < -9) {
            bpEnd = $f4["bpEnd"]
            $f4["bpEnd"] = cen[chr];
            $f4["length"] = $f4["bpEnd"] - $f4["bpStart"];
            $f4["p"]="T";
            $f4["q"]="N";
            $f4["depth"]=obs_p/exp_p-1
            $f4["depthSE"]=sqrt(obs_p)/exp_p
            print $0;

            $f4["bpStart"] = cen[chr];
            $f4["bpEnd"] = bpEnd;
            $f4["length"] = $f4["bpEnd"] - $f4["bpStart"];
            $f4["p"]="N";
            $f4["q"]="T";
            $f4["depth"]=obs_q/exp_q-1
            $f4["depthSE"]=sqrt(obs_q)/exp_q
            print $0;
            skip=1;
        }
    }
    if (skip==0) {
        print $0;
    }
}
' \
    <(zcat /mnt/project/lohdata/david/mCAs_WGS/GRCh38_supp_files/cytoBand.txt.gz  | grep acen | grep p) \
    $CALLS_FILE \
    <(cat $DEPTH_PROFILES | grep -v nan ) \
    $CALLS_FILE

set -euo pipefail

awk '
BEGIN { OFS="\t" }
NR==1 || $1!="ID" {   
    for(left=-2000; left<=2000; left+=1000) {
        for(right=-2000; right<=2000; right+=1000) {
            if(($1 OFS $2 OFS $3+left OFS $4+right) in total) {
                $3=$3+left;
                $4=$4+right;
            }
        }
    }
    region=($1 OFS $2 OFS $3 OFS $4);
    left_forward[region]+=$5;
    right_forward[region]+=$6;
    total[region]+=$7;
} 
END {
    for (region in total) {
        print region,left_forward[region],right_forward[region],total[region];
    }
}
'

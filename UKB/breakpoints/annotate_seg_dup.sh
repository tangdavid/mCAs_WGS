set -euo pipefail

BREAKPOINT_FILE=$1
BED_FILE=/mnt/project/lohdata/david/mCAs_WGS/breakpoints/ucsc.SegDup.bed


sort -k1,1 -k2,2 -u $BREAKPOINT_FILE \
    | awk '
{
    seen[$1]+=1; 
    lines[$0]=$1
} 
END {
    for (line in lines) {
        if (seen[lines[line]]==1) print line
    }
}' \
    | awk -v OFS='\t' '{
        print $2,$3-5e1,$3+5e1,$1,$3,$4; 
        print $2,$4-5e1,$4+5e1,$1,$3,$4
    }' \
    | bedtools intersect \
        -a - \
        -b $BED_FILE \
        -wao \
    | awk '{
        split($10,a,":"); 
        if(a[1]!=$1) next; 
        if($5>$8 && $5<$9 && $6>a[2] && $6<a[2]+$9-$8)  print $4,$1,$5,$6,$10,$12
    }'  | column -t

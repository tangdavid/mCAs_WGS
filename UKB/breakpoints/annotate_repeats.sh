set -euo pipefail

BREAKPOINT_FILE=$1
BED_FILE=/mnt/project/lohdata/david/mCAs_WGS/breakpoints/ucsc.RepeatMasker.bed


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
        print $2,$4-5e1,$4+5e1,$1,$3,$4;
    }' \
    | bedtools intersect \
        -a - \
        -b $BED_FILE \
        -wao \
    | awk -v OFS='\t' '{
        chrom_pos = $1":"$2+5e1;
        if(chrom_pos in repeat) {
            if (first_seen[chrom_pos] == $4)
                repeat[chrom_pos] = repeat[chrom_pos]","$10;
            else 
                recurrent[chrom_pos]+=1;
        }
        else {
            repeat[chrom_pos]=$10;
            first_seen[chrom_pos]=$5;
            recurrent[chrom_pos]=1;
        }
        lines[$4 OFS $1 OFS $5 OFS $6]=1;
    }
    END {
        for (line in lines) {
            split(line, a, OFS);
            chrom_pos_left  = a[2]":"a[3];
            chrom_pos_right = a[2]":"a[4];
            print line,repeat[chrom_pos_left],recurrent[chrom_pos_left],repeat[chrom_pos_right],recurrent[chrom_pos_right];
        }
    }' 

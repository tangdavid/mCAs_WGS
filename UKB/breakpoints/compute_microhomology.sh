set -euo pipefail

SPLIT_READS_FILE=$1
WINDOW_SIZE=50

samtools faidx \
    -r <(cat $SPLIT_READS_FILE | awk -v ws=$WINDOW_SIZE '{print $3":"$4-ws"-"$4+ws; print $7":"$8-ws"-"$8+ws}' | sort -u ) \
    -n 1000 \
    ~/tmp/GRCh38_full_analysis_set_plus_decoy_hla.fa \
    | awk -v OFS='\t' -v ws=$WINDOW_SIZE -F'\t' '
function max(x, y) {return (x>y)? x: y}
ARGIND==1 && FNR%2 == 1 {key=substr($1,2)} 
ARGIND==1 && FNR%2 == 0 {sequence[key]=$1}
ARGIND==2 {
    left_seq=sequence[$3":"$4-ws"-"$4+ws];
    right_seq=sequence[$7":"$8-ws"-"$8+ws];
    
    split($11, split_seq, ":");


    if ($6~/^[0-9]*M[0-9]*[H,S]$/ && $10~/^[0-9]*[H,S][0-9]*M$/) {
        forward_seq=substr(right_seq, ws+1, ws);
        backward_seq=substr(left_seq, 1, ws);
    } else if ($10~/^[0-9]*M[0-9]*[H,S]$/ && $6~/^[0-9]*[H,S][0-9]*M$/) {
        forward_seq=substr(left_seq, ws+1, ws);
        backward_seq=substr(right_seq, 1, ws);
    } else {
        next;
    }
    homology_seq="."
    for (homology=ws; homology >= 1; homology--) {
        if (substr(backward_seq, length(backward_seq)-(homology-1), homology) == substr(forward_seq, 1, homology)) {
            homology_seq = substr(forward_seq, 1, homology)
            break;
        }
    }
    print $0,homology,homology_seq
}
' - $SPLIT_READS_FILE

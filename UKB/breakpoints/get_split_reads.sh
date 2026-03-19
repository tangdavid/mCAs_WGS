set -euo pipefail

ID=$1
CALL=$2
BUFFER=${3:-5}

DX_PROJECT="WGS_500K"

CRAM_DIR="${DX_PROJECT}:/Bulk/DRAGEN WGS/Whole genome CRAM files (DRAGEN) [500k release]/${ID:0:2}"
CRAM_URL=`dx make_download_url "${CRAM_DIR}/${ID}_24048_0_0.dragen.cram"` 
CRAI_URL=`dx make_download_url "${CRAM_DIR}/${ID}_24048_0_0.dragen.cram.crai"`

REGION=`echo $CALL | awk -v BUFFER=$BUFFER '
BEGIN {pad=BUFFER*1e3}
{
    split($0, r, ":"); 
    split(r[2], pos, "-"); 
    print r[1]":"pos[1]-pad"-"pos[1]+pad" "r[1]":"pos[2]-pad"-"pos[2]+pad
}'`
LENGTH=`echo $CALL | awk '{split($0, r, ":"); split(r[2], pos, "-"); print pos[2]-pos[1]}'`
ENDPOINTS=`echo $REGION | tr ' ' ','`

samtools view \
    -T ${TMP_DIR}/GRCh38_full_analysis_set_plus_decoy_hla.fa ${CRAM_URL} -X ${CRAI_URL} ${REGION} \
    -F 12 \
    -q 30 \
    | awk -v sample=$ID -v call=$CALL -v ENDPOINTS=$ENDPOINTS -v OFS='\t' '
BEGIN {
    split(ENDPOINTS, endpoints, ",");
}
function abs(x) { return (x>0)? x: -x}
function max(x, y) {return (x<y) ? y : x}
function min(x, y) {return (x<y) ? x : y}
!and($2, 0x400)  && index($0, "SA:Z")>0{ 
    split("", sa);
    for (i=1; i<=NF; i++) {
        if(substr($i, 1, 4) =="SA:Z") {
            split(substr($i,6),sa,";"); 
            break;
        }
    }
    if(length(sa) > 2) next; # more than 2 SA
    split(sa[1], aln, ","); 
    if($5+aln[5]<90) next;   # combined mapping quality of at least 90
    if(!($6~/^[0-9]*M[0-9]*[H,S]$/ || $6 ~/^[0-9]*[H,S][0-9]*M$/)) next;

    pa_window=0;
    sa_windon=0;
    for (i=1; i<=length(endpoints); i++) {
        split(endpoints[i], a, ":");
        chrom=a[1];
        split(a[2],window, "-");
        if ($3==chrom && window[1] < $4 && $4 < window[2]) {
            pa_window=i;
        }
        if (aln[1]==chrom && window[1] < aln[2] && aln[2] < window[2]) {
            sa_window=i;
        }
    }
    if(pa_window + sa_window!=3) next;
    seen[$1]+=1;
    split($6, cigar, /[MHS]/);
    if($6~/^[0-9]*M[0-9]*[H,S]$/) {
        pos_adj = cigar[1];
    } else {
        pos_adj = 0;
    }

    if (seen[$1]==1) {
        left_chrom[$1]=$3;
        left_pos[$1]=$4 + pos_adj;
        left_mapq[$1]=$5;
        left_cigar[$1]=$6;
    }
    else {  
        right_chrom[$1]=$3;
        right_pos[$1]=$4 + pos_adj;
        right_mapq[$1]=$5;
        right_cigar[$1]=$6;
    }
    if(index($6, "H")) next;
    split_seq[$1] = substr($10, cigar[1]-20+1, 20)":"substr($10, cigar[1]+1, 20);
}
END {
    for (rn in seen) {
        if (seen[rn]!=2) continue; 
        print sample, call, left_chrom[rn], left_pos[rn], left_mapq[rn], left_cigar[rn], right_chrom[rn], right_pos[rn], right_mapq[rn], right_cigar[rn], split_seq[rn];
    }
}
'

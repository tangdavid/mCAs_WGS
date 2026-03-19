# set on input
# REF_DIR
#
set -euo pipefail 

ID=$1
ENDPOINTS=$2

DX_PROJECT="WGS_500K"

CRAM_DIR="${DX_PROJECT}:/Bulk/DRAGEN WGS/Whole genome CRAM files (DRAGEN) [500k release]/${ID:0:2}"
CRAM_URL=`dx make_download_url "${CRAM_DIR}/${ID}_24048_0_0.dragen.cram"` 
CRAI_URL=`dx make_download_url "${CRAM_DIR}/${ID}_24048_0_0.dragen.cram.crai"`

REGION=`echo $ENDPOINTS | awk '{
    split($0, a, ","); 
    for(i=1; i<=length(a); i++) {
        printf a[i]" ";
    }
}'
`

samtools view \
    -T ${REF_DIR}/GRCh38_full_analysis_set_plus_decoy_hla.fa ${CRAM_URL} -X ${CRAI_URL} ${REGION} \
    -F 12 \
    -q 30 \
    | cut -f-9 \
    | awk -v sample=$ID -v ENDPOINTS=$ENDPOINTS -v OFS='\t' '
function abs(x) { return (x>0)? x: -x}
BEGIN {
    split(ENDPOINTS, endpoints, ","); 
}
{ 
    second_pos="none";
    first_pos="none";
    $7=($7=="=") ? $3 : $7;
    for (i=1; i<=length(endpoints); i++) {
        split(endpoints[i], a, ":");
        chrom=a[1];
        split(a[2],window, "-");
        mid_point = int(window[1]+(window[2]-window[1])/2);
        if ($3==chrom && window[1] < $4 && $4 < window[2]) {
            first_pos=int($4/1e3)*1e3;
            first_chrom=chrom;
        }
        if ($7==chrom && window[1] < $8 && $8 < window[2]) {
            second_pos=int($8/1e3)*1e3;
            second_chrom=chrom;
        }
    }       
    
    if (first_pos=="none" || second_pos=="none") next;
    if (first_chrom==second_chrom && abs(first_pos-second_pos)< 5e3) next;
    first_pos=first_chrom OFS first_pos;
    second_pos=second_chrom OFS second_pos;
    if($1 in seen) next;
    seen[$1]=1;
    if ((second_pos OFS first_pos) in tot_support) {
        # flip the pair to merge data
        pair = second_pos OFS first_pos;
        mate_forward[pair]+=xor(rshift(and($2, 0x10), 4), 1);
        read_forward[pair]+=xor(rshift(and($2, 0x20), 5), 1);
        tot_support[pair]+=1;
    } else {
        pair=first_pos OFS second_pos;
        read_forward[pair]+=xor(rshift(and($2, 0x10), 4), 1);
        mate_forward[pair]+=xor(rshift(and($2, 0x20), 5), 1);
        tot_support[pair]+=1;
    }
}
END {
    for(pair in tot_support) {
        print sample, pair, read_forward[pair], mate_forward[pair], tot_support[pair];
    }
}
'

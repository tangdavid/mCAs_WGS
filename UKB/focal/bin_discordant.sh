# set on input:
# DISCORDANT_FILE
# REGION
# MAX_LENGTH
# MIN_LENGTH

set -euo pipefail

CHROM=`echo $REGION | cut -f1 -d":"`

grep -Fw $CHROM $DISCORDANT_FILE \
    | awk -v region=$REGION -v max_length=$MAX_LENGTH -v min_length=$MIN_LENGTH '
function print_binned_calls(seen) {
    calls = 0
    for (call in seen) {
        if(seen[call]>0 && first_forward[call]==seen[call] && second_forward[call]==0) {
            print call,seen[call],"INNER",region;
            calls += 1;
        }
    }
} 
BEGIN {
    split(region, a, ":")
    chrom=a[1]
    split(a[2], pos, "-")
    startbp=pos[1]
    endbp=pos[2]
}
prevID && $1!=prevID {
    print_binned_calls(seen); 
    delete seen
    delete first_forward;
    delete second_forward;
} 
{
    prevID=$1
}
($4>startbp && $3<endbp) && $4-$3 < max_length * 1e3 && $4-$3 > min_length * 1e3 {
    left=int($3/1e3)*1e3; 
    right=int($4/1e3)*1e3; 
    call = $1 OFS chrom":"left"-"right
    seen[call]+=1; 
    first_forward[call]+=!and($5, 0x10)
    second_forward[call]+=!and($6, 0x10)
} 
END {
    print_binned_calls(seen)
}'

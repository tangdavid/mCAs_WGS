set -euo pipefail

ANNOTATION_FILE=$1
REFERENCE_FILE=$2

awk -v OFS='\t' '
function abs(x) {return (x>0) ? x : -x} 
ARGIND==1 && FNR==1 { 
    for(i=1; i<=NF; i++) f1[$i]=i
}
ARGIND==1 {
    ID=$(f1["ID"])
    chrom=$(f1["chr"])
    start=$(f1["bpStart"])
    end=$(f1["bpEnd"])
    region_map[ID][chrom OFS start OFS end]=1
} 
ARGIND==2 && FNR==1 {
    for(i=1; i<=NF; i++) f2[$i]=i
    print
}
ARGIND==2 && FNR>1 {
    ID=$(f2["ID"])
    if(!(ID in region_map)) {next} 
    for(region in region_map[ID]) {
        split(region, a, OFS); 
        chrom=a[1]; 
        start=a[2]; 
        end=a[3]; 
        if(chrom==$(f2["chr"]) && abs($(f2["bpStart"])-start)<50e3 && abs($(f2["bpEnd"])-end)<50e3) {
            $(f2["bpStart"])=start; 
            $(f2["bpEnd"])=end; 
            print
        }
    }
}' $ANNOTATION_FILE $REFERENCE_FILE

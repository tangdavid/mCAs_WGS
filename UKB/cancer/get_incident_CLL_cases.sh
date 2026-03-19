set -euo pipefail

CANCER_PHENO_FILE=$1

cat $CANCER_PHENO_FILE \
    | awk -v OFS='\t' '
FNR==1 {
    print $0
} 
FNR>1 && $5<($1 in earliest) ? earliest[$1] : "2024" {
    earliest[$1]=$5; 
    data[$1]=$0
} 
END {
    censor_date="2020-12-31"
    for (ID in data) {
        split(data[ID], row, OFS)
    
        blood_draw = row[4]
        diagnosis = row[5]
        gsub(/-/," ",blood_draw); t1=mktime(blood_draw" 0 0 0")
        gsub(/-/," ",diagnosis); t2=mktime(diagnosis" 0 0 0")
        
        if((t2-t1)/86400/365.25 < 1 || $5>censor_date) continue
        if (row[2]==9823 && row[3]==3) # https://biobank.ndph.ox.ac.uk/showcase/ukb/docs/ICDcancermorph.pdf
            print data[ID]
    }
}'

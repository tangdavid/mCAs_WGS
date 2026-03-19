# set on input:
# CHR

REMOVE_LOCI=${REMOVE_LOCI:-sig}

set -euo pipefail
echo "`date`: removing CN-LOH explained by rare variants at ${REMOVE_LOCI} loci on $CHR..." > /dev/stderr

zcat /mnt/project/lohdata/david/mCAs_WGS/cisGWAS/results/$CHR.LoF.protein_coding.CN-LOH.GWAS.txt.gz \
    | awk '$10<1.20152e-05 && $4=="ref"' \
    | sort -k 6,6nr \
    | awk -v OFS='\t' -v LOCI=$REMOVE_LOCI '
BEGIN {
    split(LOCI, a, "|");
    for(i in a) loci[a[i]]=1;
}
{
    split($3,a,"."); 
    gene=a[1]; 
    if(gene in seen) next; 
    seen[gene]=1; 
    if (LOCI=="sig") {
        print $1,$2
    }
    else if (gene in loci) {
        print $1,$2
    }
}
END {
    print "none",1
}' \
    | bcftools query -R - /mnt/project/lohdata/resources/burden_masks/vcfs/$CHR.LoF.scaffold.phased.AS.bcf -f '[%SAMPLE\t%ID\t%GT\n]' \
    | awk -v CHR=$CHR '
ARGIND==1 && FNR==1 {for(i=1; i<=NF; i++) f[$i]=i}
ARGIND==1 && FNR>1 {
    if ($f["chr"] != CHR) next;
    if ($f["type"] != "CN-LOH") next;
    keep[$f["ID"]] = 1
}
ARGIND==2 && ($3=="0|1" || $3=="1|0") && ($1 in keep) {
    print $1
}' \
    /mnt/project/lohdata/david/mCAs_WGS/calls/WGS_500k.calls.txt \
    - \
    | sort -u 

if [[ $CHR == "chr1" && ( $REMOVE_LOCI == "sig" || $REMOVE_LOCI == "MPL" ) ]]; then
    echo "`date`: removing heterozygous carriers of chr1:43338725:G:C MPL intronic splice variant..." > /dev/stderr
    SNP_ARRAY_BCF=/mnt/project/lohdata/david/mCAs_WGS/scaffold/chr1.GRCh38.40709.AC.bcf
    bcftools query -f '[%SAMPLE\t%ID\t%GT\n]' -r chr1:43338725-43338725 $SNP_ARRAY_BCF \
        | awk '$2=="Affx-89013314" && $3!="0|0" && $3!="1|1" {print $1}' 
fi

if [[ $CHR == "chr11" && ( $REMOVE_LOCI == "sig" || $REMOVE_LOCI == "ATM" )  ]]; then 
    echo "`date`: removing heterozygous carriers of chr11:108257471:T:G ATM intronic splice variant..." > /dev/stderr
    SNP_ARRAY_BCF=/mnt/project/lohdata/david/mCAs_WGS/scaffold/chr11.GRCh38.40709.AC.bcf
    bcftools query -f '[%SAMPLE\t%ID\t%GT\n]' -r chr11:108257471-108257471 $SNP_ARRAY_BCF \
        | awk '$2=="Affx-80243477" && $3!="0|0" && $3!="1|1" {print $1}' 
fi

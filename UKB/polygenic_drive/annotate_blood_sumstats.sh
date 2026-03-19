# set on input:
# BLOOD_TRAIT

set -euo pipefail 

pip install CrossMap 2> /dev/null 1>&2
wget -nv -nc http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz 2> /dev/null

SUMSTATS_DIR=/mnt/project/lohdata/david/mCAs_WGS/polygenic_drive/sumstats
zcat $SUMSTATS_DIR/sig.${BLOOD_TRAIT}.imp_v3.stats.gz  \
    | awk '{print $2,$3,$3,$5,$6,$7,$11,$15}' \
    | CrossMap bed hg19ToHg38.over.chain.gz /dev/stdin 2> /dev/null \
    | cut -f2 -d">" \
    | awk -v OFS='\t' '
ARGIND==1 {
    # beta is flipped because BOLT puts effect allele first and I put effect allele second
    var="chr"$1 OFS $2 OFS $5 OFS $4;
    var_flipped="chr"$1 OFS $2 OFS $4 OFS $5;
    beta[var]=$7; 
    beta[var_flipped]=-$7;
    chi2[var]=$8;
    chi2[var_flipped]=$8;
    MAF[var]=($6<0.5) ? $6:1-$6;
    MAF[var_flipped]=MAF[var];
}
ARGIND==2 && FNR==1 {
    print $0,"beta","chi2","MAF";
}
ARGIND==2 && FNR>1 {
    var=$1 OFS $2 OFS $4 OFS $5
    if(var in beta) print $0,beta[var],chi2[var],MAF[var];
}
' \
    - \
    <(zcat /mnt/project/lohdata/david/mCAs_WGS/cisGWAS/results/common_var/restricted/chr{1..22}.common_var.binom_p.GWAS.txt.gz) \
    | sort -k12,12rg \
    | awk '
{
    if($1 in kept) 
        for(i in kept[$1]) 
            if($2 > i-1e6 && $2 < i+1e6) 
                next; 
    kept[$1][$2]=$0; 
    print $0;
}' \
    | sort -k1,1V -k2,2n

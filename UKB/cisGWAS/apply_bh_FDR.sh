set -euo pipefail

TMP_DIR=$HOME/tmp
mkdir -p $TMP_DIR

RESULTS_DIR=/mnt/project/lohdata/david/mCAs_WGS/cisGWAS/results
SUFFIX=LoF.protein_coding.CN-LOH.GWAS.txt.gz

zcat $RESULTS_DIR/chr{1..22}.$SUFFIX \
    | awk 'NR==1 || $4=="ref"' \
    | sort -k10,10g \
    | awk -v OFS='\t' '
BEGIN {
    split("PTPRF|MMACHC|CYP4A11|KPNA3|SPRYD7|DLEU7|KCNRG|RNASEH2B|TRIM13|OR4F15|OR4F6|TARS3",a,"|")
    for(i in a) blacklist[a[i]]=1;
    k=0;
}
NR==1 {
    for(i=1; i<=NF; i++) f[$i]=i
}
$f["ref"]=="ref" && NR>1 {
    split($f["variant"],a,".");
    gene = a[1];
    if(gene in blacklist) next;
    k++;
    pval[k]=$f["fisher_p"];
    row[k]=$f["chr"] OFS $f["pos"]/1e6 OFS gene OFS $f["alt"] OFS $f["carrier_mCA"] OFS $f["noncarrier_mCA"] OFS $f["carrier_control"] OFS $f["noncarrier_control"] OFS $f["fisher_p"];
    gene_list[k]=gene;
}
END {
    alpha = 0.01;
    m = k;
    thresh = 0;
    print "chr","pos","gene","mask","carrier_mCA","noncarrier_mCA","carrier_control","noncarrier_control","pval","threshold"
    for(k=1; k<=m; k++) {
        thresh = (pval[k] <= k/m * alpha) ? pval[k] : thresh;
    }
    for(k=1; k<=m; k++) {
        if (gene_list[k] in seen) continue;
        if (pval[k] <= thresh)  {
            seen[gene_list[k]] = 1;
            print row[k],thresh;
        }
    }
}' \
    | sort -k1,1V -k2,2n

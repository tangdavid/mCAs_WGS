set -euo pipefail 

awk -F, -v OFS='\t' '
ARGIND==1 && ($6=="Lymphoid" || $6=="Myeloid") {
    celltype[$1]=$6
} 
ARGIND==2 && FNR==1 {
    for(i=2; i<=NF; i++) {
        split($i,a," "); 
        gene[i]=a[1];
        gene_idx[a[1]]=i;
    }
} 
ARGIND==2 && FNR>1 && ($1 in celltype) {
    for(i=2; i<=NF; i++) {
        if (length($i)==0) continue;
        if(celltype[$1]=="Myeloid") {
            blood[i]+=$i;
            blood_cnt[i]++;
            myeloid[i]+=$i; 
            myeloid_cnt[i]++;
        }
        if(celltype[$1]=="Lymphoid") {
            blood[i]+=$i;
            blood_cnt[i]++;
            lymphoid[i]+=$i;
            lymphoid_cnt[i]++;
        }

    }
} 
END {
    print "gene","lymphoid_depmap","myeloid_depmap","lymphoid_myeloid_depmap"
    for(i in gene) {
        print gene[i],lymphoid[i]/lymphoid_cnt[i],myeloid[i]/myeloid_cnt[i],blood[i]/blood_cnt[i];
    }
}' \
    /mnt/project/lohdata/david/mCAs_WGS/cisGWAS/results/refined_phase/Model.DepMap.25Q2.csv \
    /mnt/project/lohdata/david/mCAs_WGS/cisGWAS/results/refined_phase/CRISPRGeneEffect.DepMap.25Q2.csv

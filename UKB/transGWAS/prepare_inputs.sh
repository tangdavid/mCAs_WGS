set -euo pipefail

plink2 \
    --bfile /mnt/project/lohdata/david/regenie/ukb_cal_chr1_22_v2_merged \
    --maf 0.01 \
    --geno 0.10 \
    --mac 100 \
    --mind 0.1 \
    --hwe 1e-15 \
    --write-snplist \
    --write-samples \
    --no-id-header \
    --out $OUT_DIR/qc_keep

cat /mnt/project/lohdata/david/mCAs_WGS/calls/WGS_500k.calls.txt  \
    | awk -v OFS='\t' '
ARGIND==1 && $2!="WGS7" && $2!="WGS8" && $2!="WGS9" {
    samples[$1]=1
}
ARGIND==2 && FNR==1 { for(i=1;i<=NF;i++) f[$i]=i }
ARGIND==2 {
    ID=$(f["ID"])
    type=$(f["type"])
    mca[ID][type]=1
    mca[ID]["ANY"]=1
}
ARGIND==2 {any_mca[$(f["ID"])]=1}
END {
    n=split("FID,IID,CN-LOH,GAIN,LOSS,ANY", cols, ",")
    for (i=1; i<=n; i++) printf "%s%s", cols[i], (i<n? OFS : ORS)
    for (sample in samples) {
        for (i=1; i<=n; i++) printf "%s%s", (i<=2 ? sample : mca[sample][cols[i]]+0), (i<n? OFS : ORS)
    }
}
' /mnt/project/lohdata/resources/WGS_depth/sample_batches/ID_40709.WGS_group.batch.txt - \
    > $OUT_DIR/mCA_pheno.txt 

# run on O2
# set on input:

set -euo pipefail

CHR=${CHR:-$1}
OUT_DIR=${OUT_DIR:-$2}

# vcf fields (via bcftools):
#      AN-AN_non_ukb >= 750K, FILTER=pass
# vep fields (via parse_gnomAD_vep.cpp):
#      2: Consequence
#      4: SYMBOL (require non-empty)
#      5: Gene [ENSG]
#      6: Feature_type (require Transcript)
#      7: Feature [ENST]
#      8: BIOTYPE (require protein_coding)
#     26: MANE_SELECT (require non-empty)
#     33: SOURCE (require Ensembl)
#     43: LoF (require HC)
#     44: LoF_filter (require empty)
#     45: LoF_flags (require empty)

CLUSTER="${HMS_CLUSTER:-RAP}"

if [[ $CLUSTER == "o2" ]]; then
    TMP_DIR=/n/scratch/users/d/dat019
    PRIMATEAI_SCORES="/n/data1/bwh/medicine/loh/david/resources/PrimateAI-3D/PrimateAI-3D.hg38.txt.gz"
else
    TMP_DIR=$HOME/tmp
    mkdir -p $TMP_DIR
    PRIMATEAI_SCORES="/mnt/project/lohdata/ronen/resources/PrimateAI-3D_scores.csv.gz"
    g++ -O3 parse_gnomAD_vep.cpp -o parse_gnomAD_vep
fi

echo "`date`: downloading gnomAD exome vcf..."
rm -f $TMP_DIR/gnomad.exomes.v4.1.sites.$CHR.vcf.bgz
wget -nv -P $TMP_DIR https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/exomes/gnomad.exomes.v4.1.sites.$CHR.vcf.bgz

echo "`date`: extracting LoF and missense variants..."
bcftools query -i 'AN-AN_non_ukb>750000' -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t%INFO/AN\t%INFO/AN_non_ukb\t%INFO/AN_non_ukb_nfe\t%INFO/AC\t%INFO/AC_non_ukb\t%INFO/AC_non_ukb_nfe\t%INFO/vep\n' \
    $TMP_DIR/gnomad.exomes.v4.1.sites.$CHR.vcf.bgz \
    | ./parse_gnomAD_vep \
    | awk -v OFS='\t' '
    ARGIND==1 && NF==14 {
        chr=$1;pos=$2;ref=$4;alt=$5; 
        if(($14 == "LoF" || $11 == "missense_variant") && $13=="MANE_SELECT") {
            missense_LoF[chr,pos,ref,alt]=$0
        }
        high_moderate_vep[chr,pos,ref,alt]=$0
    } 
    ARGIND==2 && FNR==1 {for(i=1; i<=NF; i++) f[$i]=i}
    ARGIND==2 && FNR>1 {
        chr=$f["chr"];pos=$f["pos"];ref=$f["non_flipped_ref"];alt=$f["non_flipped_alt"]; 
        if((chr,pos,ref,alt) in high_moderate_vep) {
            primateAI[chr,pos,ref,alt]=$f["score_PAI3D"];
            prediction[chr,pos,ref,alt]=$f["prediction"];
        }
    } 
    END {
        for(var in high_moderate_vep) {
            score = (var in primateAI) ? primateAI[var]:".";
            print (var in missense_LoF) ? missense_LoF[var] : high_moderate_vep[var],score;
        }
     }' \
        - \
        <(zcat $PRIMATEAI_SCORES | egrep -w "score_PAI3D|$CHR") \
    | sort -k2,2n \
    | gzip > $OUT_DIR/gnomAD.protein_coding.high_moderate_vep.$CHR.txt.gz

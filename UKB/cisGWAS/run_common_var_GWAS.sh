# set on input
# CHR
# REMOVE_LOCI

set -euo pipefail

set +e 
sudo apt-get update 2> /dev/null 1>&2
set -e
sudo apt --yes install libdeflate-dev 2> /dev/null 1>&2
sudo apt --yes install parallel 2> /dev/null 1>&2

pip3 install scipy 2> /dev/null 1>&2

TMP_DIR=$HOME/tmp
mkdir -p $TMP_DIR

cat /mnt/project/lohdata/resources/BOLT-LMM/{remove.nonEUR.FID_IID.40709.txt,w40709_CURRENT.FID_IID.txt} \
    | awk '{print $1}' \
    > $TMP_DIR/remove.samples.txt 

CHR=$CHR bash remove_related_individuals.sh >> $TMP_DIR/remove.samples.txt

if [[ -v REMOVE_LOCI ]]; then 
    REMOVE_LOCI=$REMOVE_LOCI CHR=$CHR bash remove_rare_variant_carriers.sh | sort -u >> $TMP_DIR/remove.samples.txt 
fi

echo "`date`: extracting cases and controls for $CHR..."

export CALLS_FILE=/mnt/project/lohdata/david/mCAs_WGS/calls/WGS_500k.calls.txt

extract_case_control() {
    read -r WGS_GROUP BATCH <<< "$1"
    INPUT_VCF=/mnt/project/lohdata/resources/imputedCommon/$WGS_GROUP/$WGS_GROUP.batch$BATCH.imputedCommon.bcf
    echo $INPUT_VCF $CHR 1>&2
    bcftools view -r $CHR $INPUT_VCF -Ou 2> /dev/null \
        | ../../bin/mCAs_WGS compute-cis-fisher \
            --calls $CALLS_FILE \
            --remove $TMP_DIR/remove.samples.txt \
            2> /dev/null \
            > $TMP_DIR/$WGS_GROUP.batch$BATCH.$CHR.case_control.txt 
}

export -f extract_case_control
export CHR TMP_DIR

cat /mnt/project/lohdata/resources/WGS_depth/sample_batches/WGS_group.batch.Nsamples.txt \
    | awk '$1!="WGS7" && $1!="WGS8" && $1!="WGS9" {print $1, $2}' \
    | parallel --halt-on-error 2 --joblog=/dev/stderr extract_case_control


echo "`date`: merging case control numbers..."
cat $TMP_DIR/WGS*.$CHR.case_control.txt \
    | awk -v OFS='\t' '
NR==1 {
    for(i=1; i<=NF; i++) f[$i]=i
    print "chr","pos","variant","ref","alt","carrier_mCA","noncarrier_mCA","carrier_control","noncarrier_control";
} 
NR>1 && $1!="chr" {
    variant=$f["chr"] OFS $f["pos"] OFS $f["variant"] OFS $f["ref"] OFS $f["alt"];
    c_mCA[variant]=c_mCA[variant]+$f["carrier_mCA"];
    nc_mCA[variant]=nc_mCA[variant]+$f["noncarrier_mCA"];
    c_control[variant]=c_control[variant]+$f["carrier_control"];
    nc_control[variant]=nc_control[variant]+$f["noncarrier_control"];
}
END {
    for (variant in c_mCA) {
        print variant, c_mCA[variant], nc_mCA[variant], c_control[variant], nc_control[variant];
    }
}
' \
    | sort -k1,1 -k2,2n \
    > $TMP_DIR/$CHR.case_control.txt 

echo "`date`: computing fisher p values for $CHR..."

python3 - > $OUT_DIR/$CHR.common_var.fisher_p.GWAS.txt <<EOF
from scipy.stats import fisher_exact

def fisher_test(c_mCA, nc_mCA, c_control, nc_control):
    contingency_table = [[c_mCA, nc_mCA],
                         [c_control, nc_control]]
    _, p_value = fisher_exact(contingency_table)
    return p_value

col2idx = {}

with open("$TMP_DIR/$CHR.case_control.txt") as f:
    for i, line in enumerate(f):
        line = line.strip()
        if (i == 0): 
            for j,col in enumerate(line.split("\t")):
                col2idx[col] = j
            print(line + "\t" + "fisher_p")
        else: 
            line_vec = line.split("\t")
            c_mCA = int(line_vec[col2idx["carrier_mCA"]])
            nc_mCA = int(line_vec[col2idx["noncarrier_mCA"]])
            c_control = int(line_vec[col2idx["carrier_control"]])
            nc_control = int(line_vec[col2idx["noncarrier_control"]])
            p_val = fisher_test(c_mCA, nc_mCA, c_control, nc_control)
            print(line + "\t" + str(p_val))

EOF

ALLELIC_SHIFT_VCF=/mnt/project/lohdata/david/mCAs_WGS/calls/allelic_shift/WGS_500k.$CHR.allelic_shift.bcf
../../bin/mCAs_WGS compute-cis-fisher \
    --calls $CALLS_FILE \
    --remove $TMP_DIR/remove.samples.txt \
    < $ALLELIC_SHIFT_VCF \
    | cut -f-5,11- \
    > $OUT_DIR/$CHR.common_var.binom_p.GWAS.txt

gzip $OUT_DIR/$CHR.common_var.binom_p.GWAS.txt
gzip $OUT_DIR/$CHR.common_var.fisher_p.GWAS.txt

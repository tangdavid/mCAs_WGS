# set on input
# TMP_DIR

set -euo pipefail

SAMPLE=$1
REGION=$2
MOCHA_FILE=$3
REF_BIAS_FILE=$4

CHR=`echo $REGION | cut -f1 -d':'`

bcftools view -s ${SAMPLE} -r ${REGION} ${MOCHA_FILE} \
    | bcftools view -g het -Ob \
    | tee ${TMP_DIR}/${SAMPLE}.${REGION}.bcf \
    | bcftools query -f '%POS\t%REF\t%FIRST_ALT\t[%AS]\n' \
    > ${TMP_DIR}/${SAMPLE}.${REGION}.AS.txt 
bcftools index -f ${TMP_DIR}/${SAMPLE}.${REGION}.bcf

NGS=WES TMP_DIR=${HOME}/tmp bash ../../tools/RAP/extract_cram_AD_RAP.sh $SAMPLE $REGION ${TMP_DIR}/${SAMPLE}.${REGION}.bcf \
    | awk -v OFS='\t' -v chrom=$CHR -v sample=${SAMPLE} -v region=${REGION} '
    ARGIND==1 && $5+$6>1e3 && $1==chrom {
        bias[$2 OFS $3 OFS $4] = $5/($5+$6);
    }
    ARGIND==2 {
        var = $2 OFS $3 OFS $4;
        if (var in bias) {
            p=bias[var];
            if (p > 0.6 || p < 0.4) {
                next;
            }
            if (p > 0.5) {
                ref_adjusted[var]=int($5 * (1-p)/p + 0.5);
                alt_adjusted[var]=$6;
                
                ref_variance[var]=($5+$6)*(1-p)*(1-p)*(1-p) / p;
                alt_variance[var]=($5+$6)*(1-p)*p;
            } else {
                ref_adjusted[var]=$5;
                alt_adjusted[var]=int($6 * p/(1-p) + 0.5);

                ref_variance[var]=($5+$6)*(1-p)*p;
                alt_variance[var]=($5+$6)*p*p*p/(1-p);
            }
        }
    }
    ARGIND==3 {
        var = $1 OFS $2 OFS $3;
        allelic_shift = $4;
        if (allelic_shift > 0) {
            major += ref_adjusted[var];
            minor += alt_adjusted[var];
            variance += ref_variance[var];
        } else if (allelic_shift < 0) {
            minor += ref_adjusted[var];
            major += alt_adjusted[var];
            variance += alt_variance[var];

        }
    }
    END {  
        n = major + minor;
        if (n == 0 ) {
            z = 0;
        } else {
            z = (major - 0.5 * n) / sqrt(variance);
        }
        print sample, region, major, minor, z;
    }
' \
    ${REF_BIAS_FILE} \
    - \
     ${TMP_DIR}/${SAMPLE}.${REGION}.AS.txt

rm ${TMP_DIR}/${SAMPLE}.${REGION}.bcf{,.csi}
rm ${TMP_DIR}/${SAMPLE}.${REGION}.AS.txt

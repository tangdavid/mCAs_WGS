# set on input
# TMP_DIR

set -euo pipefail

CALLS_FILE=$1
MOCHA_FILE=$2
OUT_FILE=$3
OUT_FILE_BASE=`basename $OUT_FILE`

source /mnt/project/lohdata/david/tools/python_setup_RAP.sh 2> /dev/null 1>&2

REF_FILE=GRCh38_full_analysis_set_plus_decoy_hla.fa
[[ ! -f ${TMP_DIR}/${REF_FILE} ]] && cp /mnt/project/lohdata/resources/GRCh38/${REF_FILE}{,.fai} ${TMP_DIR}

if [[ `wc -l < ${CALLS_FILE}` == 1 ]]; then
    touch ${OUT_FILE}
    exit 0
fi

cat /mnt/project/lohdata/david/mCAs_WGS/cramADs_WES/10k.WES.ref_bias.txt > ${TMP_DIR}/10k.WES.ref_bias.txt

set +e
awk '
    FNR==1 { for (i=1; i<=NF; i++) f[$i]=i; }
    FNR > 1 {print $(f["ID"])" "$(f["chr"])":"$(f["bpStart"])"-"$(f["bpEnd"])}
    ' ${CALLS_FILE} \
    | parallel --colsep ' ' --joblog /dev/stderr --retries 4 \
         TMP_DIR=$TMP_DIR bash compute_WES_AI.sh {1} {2} ${MOCHA_FILE} ${TMP_DIR}/10k.WES.ref_bias.txt > ${TMP_DIR}/${OUT_FILE_BASE}.tmp
set -e

awk -v OFS='\t' '
    ARGIND==1 && FNR==1 { for (i=1; i<=NF; i++) f[$i]=i; }
    ARGIND==1 && FNR>1 {
        call = $(f["ID"])" "$(f["chr"])":"$(f["bpStart"])"-"$(f["bpEnd"]);
        bdev[call] = $(f["bdev"]);
        bdev_se[call] = $(f["bdevSE"]);
    }
    ARGIND==2 {
        call = $1" "$2;
        print $0,bdev[call],bdev_se[call];
    }
' \
    ${CALLS_FILE} \
    ${TMP_DIR}/${OUT_FILE_BASE}.tmp \
    | python3 compute_expected_replication.py \
    | tr ':' '\t' \
    | sed 's/-/\t/' \
    > ${OUT_FILE}



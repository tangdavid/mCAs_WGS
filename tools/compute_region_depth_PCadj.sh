# set on input:
# PC_DIR
# DEPTH_DIR
# PFX 
# ID_REGION_LIST

set -euo pipefail

: ${PC_DIR:?PC_DIR is not set}
: ${DEPTH_DIR:?DEPTH_DIR is not set}
: ${PFX:?PFX is not set}
: ${ID_REGION_LIST?ID_REGION_LIST is not set}
: ${GC_MASKS?GC_MASKS is not set}


PROFILE_OUT=${1:-/dev/null}

if [[ -z "$TMP_DIR" ]]; then
    TMP_DIR=$(mktemp -d /tmp/compute_region_depth_PCadj.XXXXXX)
else
    mkdir -p "$TMP_DIR"
fi

# extract bin counts per sample

TMP_PFX=$TMP_DIR/`basename $PFX`

rm -rf $TMP_PFX.counts
unzip $DEPTH_DIR/$PFX.counts.zip -d $TMP_PFX.counts >&2 

RES_KB=${RES_KB:-20}
MAX_DIFF=${MAX_DIFF:-0.08}
NUM_PCS=${NUM_PCS:-20}
CNV_FILE=${CNV_FILE:-$DEPTH_DIR/$PFX.CNVs.txt.gz}

PC_PFX=${PC_PFX:-profiles_${RES_KB}kb.maxDiff$MAX_DIFF}

# compute GC+PC-adjusted depths in each query region
`git rev-parse --show-toplevel`/bin/computeRegionDepthsPCadj \
    $ID_REGION_LIST \
    $CNV_FILE \
    $DEPTH_DIR/$PFX.${RES_KB}kb_profiles.txt.gz \
    $PC_DIR/$PC_PFX \
    $GC_MASKS \
    $TMP_PFX.counts/%s.counts.bin.bgz \
    $MAX_DIFF \
    $NUM_PCS \
    $OUT_FILE \
    $PROFILE_OUT 

rm -rf $TMP_PFX.counts

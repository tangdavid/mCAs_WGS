# set on input: WGS_GROUP, BATCH, ID_REGION_LIST, OUT_FILE

set -euo pipefail

: ${ID_REGION_LIST?ID_REGION_LIST is not set}
: ${WGS_GROUP?WGS_GROUP is not set}
: ${BATCH?BATCH is not set}
: ${OUT_FILE?OUT_FILE is not set}

DEPTH_DIR=/mnt/project/lohdata/resources/WGS_depth
PC_DIR=$DEPTH_DIR/PCadj/WGS1-6
PFX=$WGS_GROUP/$WGS_GROUP.batch$BATCH


PROFILE_OUT=${1:-/dev/null}

TMP_DIR=$HOME/tmp
mkdir -p $TMP_DIR

# extract bin counts per sample

TMP_PFX=$TMP_DIR/`basename $PFX`

rm -rf $TMP_PFX.counts
for I in {1..10}; do # unsquashfs from /mnt/project usually needs to be retried a few times
    unsquashfs -d $TMP_PFX.counts $DEPTH_DIR/$PFX.counts.sqfs 1>&2 && break
done

RES_KB=20
MAX_DIFF=0.08
NUM_PCS=${NUM_PCS:-20}

# compute GC+PC-adjusted depths in each query region
`git rev-parse --show-toplevel`/bin/computeRegionDepthsPCadj \
    $ID_REGION_LIST \
    $DEPTH_DIR/$PFX.CNVs.txt.gz \
    $DEPTH_DIR/$PFX.${RES_KB}kb_profiles.txt.gz \
    $PC_DIR/profiles_${RES_KB}kb.maxDiff$MAX_DIFF \
    $DEPTH_DIR/GRCh38.GCmasks.bin \
    $TMP_PFX.counts/%s.counts.bin.bgz \
    $MAX_DIFF \
    $NUM_PCS \
    $OUT_FILE \
    $PROFILE_OUT 

rm -rf $TMP_PFX.counts

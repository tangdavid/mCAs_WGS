# set on input: $WGS_GROUP $BATCH (smallest: WGS2 58 = 115 samples)

set -euo pipefail

OUT_PREFIX=$WGS_GROUP.batch$BATCH

set +e # OK if apt-get update fails
sudo apt-get update
set -e
sudo apt --yes install parallel
sudo apt --yes install pigz

unset DX_WORKSPACE_ID # to allow dx make_download_url

TMP_DIR=$HOME/tmp
mkdir -p $TMP_DIR

# create sample list
awk -v wgs=$WGS_GROUP -v batch=$BATCH '$2==wgs && $3==batch {print $1}' \
    /mnt/project/lohdata/resources/WGS_depth/sample_batches/ID_40709.WGS_group.batch.txt \
    > $TMP_DIR/IDs.txt

wc -l $TMP_DIR/IDs.txt

run_count_reads() {
    ID=$1
    TMP_DIR=$2

    ../../../bin/countReadsSaveDiscordant \
        $ID \
        <(curl -o - `dx make_download_url "WGS_500K:/Bulk/DRAGEN WGS/Whole genome CRAM files (DRAGEN) [500k release]/${ID:0:2}/${ID}_24048_0_0.dragen.cram"` 2> /dev/null ) \
        GRCh38.GCmasks.bin \
        1000 \
        $TMP_DIR/$ID \
        > /dev/null
}
export -f run_count_reads

# count reads in each cram file
set +e # OK if some jobs fail; status will appear in job log
cat $TMP_DIR/IDs.txt \
    | parallel --joblog /dev/stderr --retries 4 run_count_reads {} $TMP_DIR
set -e

# collate GC profiles
cat $TMP_DIR/*.GCprofile.txt > $OUT_PREFIX.GCprofiles.txt

# collate 1kb-bin read counts into archive
mksquashfs $TMP_DIR/*.counts.bin.bgz $OUT_PREFIX.counts.sqfs -noDataCompression -no-fragments

# collate discordant reads
cat $TMP_DIR/*.discordant.txt | pigz > $OUT_PREFIX.discordant.txt.gz

# compute 20kb and 50kb read-depth profiles (and CNV calls) simultaneously
for RES_KB in 20 50
do
    ../../../bin/computeDepthProfiles \
	<(cat $TMP_DIR/*.counts.bin.bgz) \
	GRCh38.GCmasks.bin \
	${RES_KB}e3 \
	$TMP_DIR/${RES_KB}kb && pigz $TMP_DIR/${RES_KB}kb.profiles.txt &
done
wait
for RES_KB in 20 50
do
    mv $TMP_DIR/${RES_KB}kb.profiles.txt.gz $OUT_PREFIX.${RES_KB}kb_profiles.txt.gz
    if [ ! -f $OUT_PREFIX.CNVs.txt.gz ]; then
	pigz -c $TMP_DIR/${RES_KB}kb.CNVs.txt > $OUT_PREFIX.CNVs.txt.gz
    fi
done

wc -l $OUT_PREFIX.GCprofiles.txt \
    | awk '{printf("Counted reads and created depth profiles for %d ",$1)}'
wc -l $TMP_DIR/IDs.txt \
    | awk '{printf(" of %d cram files\n",$1)}'

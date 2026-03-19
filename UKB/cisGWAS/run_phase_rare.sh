# set on input: 
# OUT_DIR
# CHR
# BATCH

set -euo pipefail

set +e 
sudo apt-get update
set -e
sudo apt --yes install parallel


GENETIC_MAP="/mnt/project/lohdata/david/mCAs_WGS/GRCh38_supp_files/genetic_map.$CHR.txt"
CHUNKS="/mnt/project/lohdata/david/mCAs_WGS/imputation_reference/small_chunks_6cM/jobs/chunks_$CHR.$BATCH.txt"

TMP_DIR=$HOME/tmp
mkdir -p $TMP_DIR

echo "`date`: downloading input bcfs..."
rm -f $TMP_DIR/$CHR.{phase_rare.input,snp_scaffold.subset_IDs}.bcf{,.csi}
cp /mnt/project/lohdata/resources/burden_masks/vcfs/tmp/$CHR.{phase_rare.input,snp_scaffold.subset_IDs}.bcf{,.csi} $TMP_DIR/

wget -nv https://github.com/odelaneau/shapeit5/releases/download/v5.1.1/phase_rare_static -O phase_rare_static && chmod u+x phase_rare_static

# Create a function to process each chunk
process_chunk() {
    read -r CHUNK CHR_NAME SCAFFOLD_REGION INPUT_REGION REMAINING <<< "$1"
    echo "`date`: phasing chunk $CHR_NAME.$CHUNK: $INPUT_REGION..."
    if ! bcftools view $TMP_DIR/$CHR.phase_rare.input.bcf -r $INPUT_REGION | grep -q '0/1'; then  
        echo "No unphased variants found; skipping chunk $CHR.$CHUNK..."
    else 
        ./phase_rare_static \
            --input $TMP_DIR/$CHR.phase_rare.input.bcf \
            --scaffold $TMP_DIR/$CHR.snp_scaffold.subset_IDs.bcf \
            --input-region $INPUT_REGION \
            --scaffold-region $SCAFFOLD_REGION \
            --output $TMP_DIR/$CHR.LoF.scaffold.phased.$CHUNK.bcf \
            --map $GENETIC_MAP \
            --thread 4 
    fi
}

export -f process_chunk
export TMP_DIR CHR GENETIC_MAP CONCAT_FILE

JOBS=$(($(nproc)/4))
if [[ $JOBS -lt 1 ]]; then JOBS=1; fi
echo $JOBS

# Run chunks in parallel
cat $CHUNKS | parallel --halt-on-error 2 --joblog /dev/stderr -j $JOBS process_chunk 

CONCAT_FILE=$TMP_DIR/$CHR.phase_rare.concat.txt
rm -f $CONCAT_FILE

while read CHUNK_LINE; do
    read -r CHUNK REMAINING <<< "$CHUNK_LINE"
    if [[ -f $TMP_DIR/$CHR.LoF.scaffold.phased.$CHUNK.bcf ]]; then
        echo $TMP_DIR/$CHR.LoF.scaffold.phased.$CHUNK.bcf >> $CONCAT_FILE
    fi
done < $CHUNKS 

mkdir -p $OUT_DIR/tmp
echo "`date`: concatenating phased chunks..."
if [[ -f $CONCAT_FILE ]]; then
    bcftools concat -n -f $CONCAT_FILE -Ob -o $OUT_DIR/tmp/$CHR.LoF.scaffold.phased.b${BATCH}.bcf
fi 

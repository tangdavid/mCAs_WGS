set -euo pipefail

FIMO_OUT=~/out_dir/fimo_output
mkdir -p ~/out_dir/fimo_output
chmod 777 $FIMO_OUT

docker run \
    -v ~/out_dir/RAG-motif_hepnon_combined.meme.txt:/data/RAG-motif_hepnon_combined.meme.txt \
    -v ~/tmp/sequence_context.fasta:/data/sequence_context.fasta \
    -v $FIMO_OUT:/output  \
    memesuite/memesuite fimo --oc /output --thresh 1.0E-4 /data/RAG-motif_hepnon_combined.meme.txt /data/sequence_context.fasta

chmod 775 $FIMO_OUT


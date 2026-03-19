set -euo pipefail 

pip3 install numpy scipy 2> /dev/null 1>&2

zcat /mnt/project/lohdata/david/resources/chr21.afreq.txt.gz | tail +2 | shuf -n 100000 | gzip > $TMP_DIR/chr21.afreq.txt.gz

awk -v OFS='\t' '
BEGIN {
    print -0.2,0.25,100,0.002
    print -0.2,0.5,100,0.002
    print -0.2,1,100,0.002

    print -0.5,0.25,100,0.002
    print -0.5,0.5,100,0.002
    print -0.5,1,100,0.002

    print -0.2,0.25,100,0.02
    print -0.2,0.5,100,0.02
    print -0.2,1,100,0.02

    print -0.5,0.25,100,0.02
    print -0.5,0.5,100,0.02
    print -0.5,1,100,0.02

    print -0.2,0.25,1000,0.002
    print -0.2,0.5,1000,0.002
    print -0.2,1,1000,0.002

    print -0.5,0.25,1000,0.002
    print -0.5,0.5,1000,0.002
    print -0.5,1,1000,0.002
}' \
    | while read LINE; do
    read -r ALPHA H2 NUM_CAUSAL PREVALENCE <<< $LINE
    python3 simulate_genetic_architecture.py $ALPHA $H2 $NUM_CAUSAL $PREVALENCE --n 50000 --AF $TMP_DIR/chr21.afreq.txt.gz
done \
    | awk 'NR==1 || $1!="ALPHA"' \
    > $OUT_DIR/genetic_architecture_simulations.txt

set -euo pipefail

export TMP_DIR=$HOME/tmp
mkdir -p $TMP_DIR

cp /mnt/project/lohdata/david/resources/gcta64 ./ && chmod 744 gcta64


echo -e "reml\treml_se\th2\tn\tm" > $OUT_FILE

run_simulations() {
    python3 simulate_h2_estimation.py --h2 $H2 --N $N --M $M --seed $SEED --out-dir $TMP_DIR
    PFX=$TMP_DIR/${N}_${M}_${H2}_${SEED}
    ./gcta64 --reml --pheno $PFX.pheno --grm $PFX --prevalence 0.5 --reml-no-constrain \
        | grep 'V(G)/Vp_L' | awk -v n=$N -v m=$M -v h2=$H2 -v OFS='\t' '{print $2,$3,h2,n,m}'
    rm $PFX.{grm.bin,grm.id,pheno}
}

export -f run_simulations

#H2=0.1
#N=500
#M=5000

set +e
echo {0.1,0.3,0.5,0.7,0.9},{100,500,1000},{1000,5000,10000},{1..100} \
    | tr ' ' '\n' \
    | parallel --colsep=, --joblog=/dev/stderr "H2={1} N={2} M={3} SEED={4} run_simulations" \
    >> $OUT_FILE
set -e



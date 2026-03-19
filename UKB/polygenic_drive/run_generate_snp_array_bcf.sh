# set on input
# OUT_DIR

set -euo pipefail

set +e 
sudo apt-get update
set -e
sudo apt --yes install parallel

echo chr{1..22} \
    | tr ' ' '\n' \
    | parallel --joblog=/dev/stderr OUT_DIR=$OUT_DIR bash generate_snp_array_bcf.sh {}

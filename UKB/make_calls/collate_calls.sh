set -euo pipefail


TYPE=$1
OUT_FILE=$2

zcat /mnt/project/lohdata/david/mCAs_WGS/calls/*/*.${TYPE}.txt.gz | awk 'NR==1 || $1!="ID"' | gzip > ${OUT_FILE}

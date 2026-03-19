set -euo pipefail

cat /mnt/project/lohdata/david/mCAs_WGS/calls/WGS_500k.calls.txt \
    | awk '$7!="T" && $8!="T"' \
    | awk '{print $1,$3,$4,$5}' \
    | bash extract_discordant_reads.sh /dev/stdin \
    | gzip > $OUT_DIR/WGS_500k.calls.all_discordant_reads.txt.gz

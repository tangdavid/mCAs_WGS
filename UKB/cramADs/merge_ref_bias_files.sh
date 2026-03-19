set -euo pipefail
NGS=WGS
DIR=/mnt/project/lohdata/david/mCAs_WGS/cramADs_${NGS}

for CHR in chr{1..22} chrX; do
    echo ${CHR} > /dev/stderr 
    cat ${DIR}/*/*.${NGS}.ref_bias.txt \
        | sed 's/chr23/chrX/g' \
        | grep ${CHR} \
        | awk -v OFS='\t' -v CHR=${CHR} '
            $1==CHR {
                var=$2 OFS $3 OFS $4; 
                ref[var] += $5; 
                alt[var] += $6
            } 
            END { 
                for (key in ref) { 
                    print CHR OFS key OFS ref[key] OFS alt[key] 
                } 
            }' \
        | sort -k2,2n
done

# set on input: 
# TMP_DIR

set -euo pipefail

if [ "$#" -lt 5 ]; then 
    echo "Usage: bash $0 ID REGION BCF_FILE REF CRAM_URL [CRAI_URL]"
    exit 1
fi

ID=$1
REGION=$2
BCF_FILE=$3
REF=$4
CRAM_URL=$5
CRAI_URL=${6:-$CRAM_URL.crai}

TMP_DIR=${TMP_DIR:-/tmp}

CHR=`echo $REGION | cut -f1 -d':'`

PFX=${TMP_DIR}/${ID}.${REGION}

# get the phased genotypes and mask out QC'd sites
bcftools view ${BCF_FILE} -s ${ID} -r ${REGION} -Ou \
    | bcftools view --no-version -Ov -g 'het' -p \
    > ${PFX}.hets.vcf

# compute mpileup; provide as 2nd input to subsequent awk post-processing command
samtools mpileup \
    -l <(awk '!index($1,"#") {print $1,$2}' ${PFX}.hets.vcf) \
    -r $REGION \
    -f $REF \
    --output-QNAME \
    --no-BAQ \
    ${CRAM_URL} -X ${CRAI_URL} \
    | awk '$6!="*"' \
    | cut -f2,5,7 \
    | sed "s/\^.//g" \
    | tr -d '$' \
    | perl -lane '$F[1] =~ s/(.([+-](\d+)(??{".{$3}"}))*)/ $1/g; print "@F"' \
    | awk -v OFS=$'\t' -v CHR=$CHR '
ARGIND==1 {
  pos2ref_alt[$2]=$4"_"$5;   # set up map from het positions to REF_ALT alleles
}
ARGIND==2 {
  POS = $1;
  if (POS - prevPOS > 5000) delete QNAME_used;
  prevPOS = POS;
  split(pos2ref_alt[POS],alleles,"_");  # look up REF and ALT for het call at POS

  REF = alleles[1];
  REF_mpileup_1 = ".";
  REF_mpileup_2 = ",";

  ALT = alleles[2];
  if (length(REF)==1 && length(ALT)==1) { # variant is SNP
    ALT_mpileup_1 = ALT;
    ALT_mpileup_2 = tolower(ALT);
  }
  else if (length(REF)==1) {              # variant is INS
    ALT_mpileup_1 = ".+"(length(ALT)-1)substr(ALT,2);
    ALT_mpileup_2 = ",+"(length(ALT)-1)tolower(substr(ALT,2));
  }
  else if (length(ALT)==1) {              # variant is DEL
    ALT_mpileup_1 = ".-"(length(REF)-1)substr(REF,2);
    ALT_mpileup_2 = ",-"(length(REF)-1)tolower(substr(REF,2));
  }
  else
    next;

  split($NF,QNAMEs,",");
  if (length(QNAMEs) != NF-2) {
    print "Error parsing mpileup: ",$0 > "/dev/stderr";
    next;
  }

  REF_AD_uniq = 0; ALT_AD_uniq = 0; REF_AD = 0; ALT_AD = 0; 
  for (i=1; i<=NF-2; i++) {
    bases = $(i+1);
    is_REF = (bases==REF_mpileup_1 || bases==REF_mpileup_2);
    is_ALT = (bases==ALT_mpileup_1 || bases==ALT_mpileup_2);
    if (is_REF != is_ALT) { # fragment supports either REF or ALT => increment AD
      if (is_REF) REF_AD++;
      else ALT_AD++;
      if (!QNAME_used[QNAMEs[i]]) { # fragment has not yet been seen => increment AD_uniq
        QNAME_used[QNAMEs[i]] = 1;
        if (is_REF) REF_AD_uniq++;
        else ALT_AD_uniq++;
      }
    }
  }

  print CHR,POS,REF,ALT,REF_AD_uniq,ALT_AD_uniq,REF_AD,ALT_AD;
}' \
    <(bcftools view -H ${PFX}.hets.vcf) \
    - 

rm ${PFX}.hets.vcf

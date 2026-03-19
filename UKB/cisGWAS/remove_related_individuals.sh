# set on input:
# CHR

set -euo pipefail

BCF_FILE=${BCF_FILE:-/mnt/project/lohdata/resources/burden_masks/vcfs/$CHR.LoF.scaffold.phased.AS.bcf}
CALL_FILE=${CALL_FILE:-/mnt/project/lohdata/david/mCAs_WGS/calls/WGS_500k.calls.txt}
CN=${CN:-"CN-LOH"}

awk -v CHR=$CHR -v CN=$CN '
ARGIND==1 {analyzed[$1]=1; }
ARGIND==2 { age[$1]=$2; }
ARGIND==3 && FNR==1 {for(i=1; i<=NF;i++) f[$i]=i}
ARGIND==3 && FNR>1 {
    if ($f["chr"] != CHR) next;
    if ($f["type"] != CN) next;
    call[$f["ID"]] = 1
}
ARGIND==4 {
    ID1 = $1;
    ID2 = $2;
    if(!(ID1 in analyzed) || !(ID2 in analyzed)) next;
    if(ID1 in remove || ID2 in remove) next;
    if (ID1 in call && ID2 in call) {
        remove[ID1] = 1;
        case_case_remove++;
    } else if (ID1 in call && !(ID2 in call)) {
        remove[ID2] = 1;
        case_ctrl_remove++;
    } else if (!(ID1 in call) && ID2 in call) {
        remove[ID1] = 1;
        case_ctrl_remove++;
    } else {
        if(age[ID1] < age[ID2]) {
            remove[ID1] = 1;
        } else {
            remove[ID2] = 1;
        }
        ctrl_ctrl_remove++;
    }
}
END {
    for(ID in remove) print ID;
    for(ID in call) {
        tot_cases++;
        removed_cases+=(ID in remove) ? 1 : 0;
    }
    print "Total cases (pre remove relateds):", tot_cases > "/dev/stderr";
    print "Total cases (post remove relateds):", tot_cases-removed_cases > "/dev/stderr";
    print "Removed (case-case):", case_case_remove+0 >> "/dev/stderr";
    print "Removed (case-ctrl):", case_ctrl_remove+0 >> "/dev/stderr";
    print "Removed (ctrl-ctrl):", ctrl_ctrl_remove+0 >> "/dev/stderr";
}
' \
    <( bcftools query -l $BCF_FILE) \
    "/mnt/project/lohdata/david/mCAs_WGS/sample_data/ID_age.40709.txt" \
    $CALL_FILE \
    "/mnt/project/Bulk/Genotype Results/Genotype calls/ukb_rel.dat"



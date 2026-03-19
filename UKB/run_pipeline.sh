set -euo pipefail

PIPELINE_DIR=$1
PIPELINE_SCRIPT=$2

mkdir -p ${HOME}/.ssh
cp /mnt/project/lohdata/david/ssh/* ${HOME}/.ssh
chmod 400 ${HOME}/.ssh/id_rsa
chmod 644 ${HOME}/.ssh/known_hosts
ssh-agent bash -c 'ssh-add ${HOME}/.ssh/id_rsa; GIT_SSH_COMMAND="ssh -o UserKnownHostsFile=${HOME}/.ssh/known_hosts" git clone -b dev git@github.com:tangdavid/mCAs_WGS.git ${HOME}/mCAs_WGS'

OUT_DIR=`pwd`
cd ${HOME}/mCAs_WGS/UKB/${PIPELINE_DIR}
echo "current commit: `git log -1 --format=%H`"
OUT_DIR=$OUT_DIR bash ${PIPELINE_SCRIPT}


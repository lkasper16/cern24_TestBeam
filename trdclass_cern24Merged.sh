#!/bin/bash


#source setup_env.sh
source /gapps/root/Linux_RHEL7-x86_64-gcc4.8.2/root-6.18.00/bin/thisroot.csh
# only root needed

RUNNUM=${1-none}
NENTRIES=${2-0}
NTREES=${3-0}
MAXEVT=${4-0}
FRSTEVT=${5-0}

if [[ ${RUNNUM} == "none" ]] ; then
    echo "================================="
    echo " Usage: ./$0 <RunNum> [nEntries] [nTrees] [Max_Events] [First_Event]"
    echo "================================="
    exit 0;
fi

echo "====>  Process MERGED RUN=$RUNNUM <=========="

root --web=off -l <<EOC
.L trdclass_cern24.C+g
trdclass_cern24 t(${RUNNUM},${NENTRIES},${NTREES},${MAXEVT},${FRSTEVT})
t.Loop()
EOC

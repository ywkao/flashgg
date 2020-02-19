#!/bin/bash
OUTPUTDIR=$1
OUTPUTFILENAME=$2
INPUTFILENAMES=$3
INDEX=$4
CMSSW_VER=$5

ARGS=$7

WD=$PWD
echo
echo
echo
rm package.tar.gz
xrdcp root://redirector.t2.ucsd.edu///store/user/smay/FCNC/MY_DIR/package.tar.gz .
export SCRAM_ARCH=slc7_amd64_gcc700
source /cvmfs/cms.cern.ch/cmsset_default.sh
tar xvf package.tar.gz
rm package.tar.gz
cp config.json $CMSSW_VER/src/flashgg/Systematics/test/MY_DIR//
cd $CMSSW_VER/src/flashgg
scramv1 b ProjectRename
scram b -j1
eval $(scram runtime -sh)
pip install --user htcondor
cd $WD
mkdir MY_DIR/
echo "ls $X509_USER_PROXY"
ls $X509_USER_PROXY
mkdir .dasmaps
mv das_maps_dbs_prod.js .dasmaps/

echo "process.source = cms.Source('PoolSource', fileNames=cms.untracked.vstring('file:$INPUTFILENAMES'))" >> "$CMSSW_BASE/src/flashgg/Systematics/test/workspaceStd.py"
COMMAND

cd MY_DIR
echo
echo
echo "Job finished with exit code ${retval}"
echo "Files in ouput folder"
ls -ltr
if [[ $retval == 0 ]]; then
    errors=""
    for file in $(find -name '*.root' -or -name '*.xml'); do
        env -i X509_USER_PROXY=${X509_USER_PROXY} gfal-copy -p -f -t 4200 --verbose $file gsiftp://gftp.t2.ucsd.edu/hadoop/cms/store/user/smay/FCNC/MY_DIR//$file
        if [[ $? != 0 ]]; then
            errors="$errors $file($?)"
        fi
    done
    if [[ -n "$errors" ]]; then
       echo "Errors while staging files"
       echo "$errors"
       exit -2
    fi
fi

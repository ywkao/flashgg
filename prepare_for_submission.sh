#!/bin/bash

file_source="/afs/cern.ch/work/y/ykao/public/VLQ_ntuple_production"

#----------------------------------------------------------------------------------------------------
# set up for reaching xml files, json, and scripts
#----------------------------------------------------------------------------------------------------
ln -s ${file_source}/json Systematics/test/json
ln -s ${file_source}/script Systematics/test/script
ln -s ${file_source}/xmlfiles Taggers/mydata

cp -p /afs/cern.ch/work/y/ykao/public/VLQ_ntuple_production/exe_for_submission/*sh Systematics/test/

ls -lhrt Systematics/test
ls -lhrt Taggers

#----------------------------------------------------------------------------------------------------
# VO set up
#----------------------------------------------------------------------------------------------------
#echo ">>> voms-proxy-init -rfc -voms cms --valid 168:00"; voms-proxy-init -rfc -voms cms --valid 168:00;
#echo ">>> voms-proxy-info"; voms-proxy-info;
#echo ""
#echo ">>> cp /tmp/x509up_u75423 ~/"; cp /tmp/x509up_u75423 ~/; ls -lhrt ~/x509up_u75423;
#echo ""
#echo ">>> export X509_USER_PROXY=~/x509up_u75423"; export X509_USER_PROXY=~/x509up_u75423; echo ${X509_USER_PROXY};

#----------------------------------------------------------------------------------------------------
# Ready for condor submission
#----------------------------------------------------------------------------------------------------
echo ">>> Ready for condor job submission"
echo ">>> type command: cd Systematics/test/"

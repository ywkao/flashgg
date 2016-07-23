# NB this command is specific to the configuration at IC and is not gaurenteed elsewhere
#LM=/afs/cern.ch/work/s/sethzenz/fromscratch107/CMSSW_8_0_8_patch1/src/flashgg/Systematics/test/Cert_271036-275783_13TeV_PromptReco_Collisions16_JSON.txt #6.26/fb
LM=/afs/cern.ch/work/s/sethzenz/fromscratch107/CMSSW_8_0_8_patch1/src/flashgg/MetaData/work/jsons/Cert_271036-276384_13TeV_PromptReco_Collisions16_JSON_NoL1T.txt #9.17/fb
fggRunJobs.py --load data_jobs.json -d data_jobs_733 -x cmsRun workspaceStd.py maxEvents=-1 -n 500 -q 1nh -D -P useAAA=0 atIC=1 doFiducial=False lumiMask=${LM}
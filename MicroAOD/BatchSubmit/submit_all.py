import sys, os
import time
import itertools
import numpy
import json

from metis.Sample import DBSSample
from metis.CondorTask import CondorTask
from metis.StatsParser import StatsParser

job_dir = "MicroAOD/RunII"
job_tag = "RunII"
exec_path = "condor_exe.sh"
tar_path = "package.tar.xz"
hadoop_path = "ttH/{0}".format(job_dir)

cmssw_ver = "CMSSW_10_5_0"

DOSKIM = True

def getArgs(pid, dataset, json):
    datasetName = dataset.split("/")[1]
    if "EGamma" in datasetName:
      datasetName += "_2018"
    args = "processType={0} datasetName={1} conditionsJSON={2}".format(pid, datasetName, json)
    args = args.replace("processType=sig_vh", "processType=AUTO")
    return args

dsdefs = []

conds_dict = {
	"datasets_RunIISummer16.json" : "MetaData/data/MetaConditions/Era2016_RR-17Jul2018_v1.json",
	"datasets_RunIIFall17.json" : "MetaData/data/MetaConditions/Era2017_RR-31Mar2018_v1.json",
	"datasets_RunIIAutumn18.json" : "MetaData/data/MetaConditions/Era2018_RR-17Sep2018_v1.json" 
}

job_jsons = ["datasets_RunIISummer16.json", "datasets_RunIIFall17.json", "datasets_RunIIAutumn18.json"]
for js in job_jsons:
    jobs = json.load(open(js))
    for pid in jobs["processes"]:
        fpo = jobs["processes"][pid]["filesPerOutput"]
        for ds in jobs["processes"][pid]["datasets"]:
            args = getArgs(pid, ds, conds_dict[js])                
	    print args
            dsdefs.append((ds, fpo, args))

total_summary = {}
while True:
    allcomplete = True
    for ds,fpo,args in dsdefs[:]:
	if "ttHJetToGG_M105_13TeV_amcatnloFXFX_madspin_pythia8" in ds:
	  print "Skipping ttH M105 because of xrootd errors"
	  continue
        sample = DBSSample( dataset=ds )
        task = CondorTask(
                sample = sample,
                open_dataset = False,
                files_per_output = fpo,
                output_name = "test_skim.root" if DOSKIM else "myMicroAODOutputFile.root",
                tag = job_tag,
    		cmssw_version = cmssw_ver,
                executable = exec_path,
                tarfile = tar_path,
                condor_submit_params = {"sites" : "T2_US_UCSD,T2_US_CALTECH,T2_US_MIT,T2_US_WISCONSIN,T2_US_Nebraska,T2_US_Purdue,T2_US_Vanderbilt,T2_US_Florida"},
                special_dir = hadoop_path,
                arguments = args.replace(" ","|")
                )
        task.process()
        allcomplete = allcomplete and task.complete()
        # save some information for the dashboard
        total_summary[ds] = task.get_task_summary()
    # parse the total summary and write out the dashboard
    StatsParser(data=total_summary, webdir="~/public_html/dump/metis_microaod_runII/").do()
    os.system("chmod -R 755 ~/public_html/dump/metis_microaod_runII")
    if allcomplete:
        print ""
        print "Job={} finished".format(job_tag)
        print ""
        break
    #print "Sleeping 300 seconds ..."
    #time.sleep(300)
    print "Sleeping 10000 seconds ..."
    time.sleep(10000)

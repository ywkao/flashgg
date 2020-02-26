import sys, os
import json
import datetime
import glob
import copy

d = datetime.datetime.today()

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--tag", help = "tag to identify this set of babies", type=str)
parser.add_argument("--years", help = "which year(s) to submit jobs for", type=str, default="2016,2017,2018")
parser.add_argument("--date", help = "date to distinguish this submission", type=str, default=d.strftime("%d-%B-%Y"))
parser.add_argument("--procs", help = "csv list of types of processes to submit", type=str, default = "sig,data,zh,fcnc")
parser.add_argument("--couplings", help = "csv list of couplings", type=str, default = "Hct,Hut")
parser.add_argument("--dry_run", help = "don't actually make tarball or submit", action="store_true")
parser.add_argument("--fcnc_wildcard", help = "wildcard to grab all fcnc samples", default = "/hadoop/cms/store/user/smay/ttH/MicroAOD/RunII/*FCNC*")
parser.add_argument("--fcnc_only", help = "only run on fcnc locally", action="store_true")
parser.add_argument("--save_tar", help = "don't remake tarball", action="store_true")
parser.add_argument("--no_syst", help = "don't run systematics", action="store_true")
args = parser.parse_args()

import parallel_utils
from flashgg_utils import *

couplings = {}
for coupling in args.couplings.split(","):
    couplings[coupling] = {}

years = args.years.split(",")
syst = "False" if args.no_syst else "True"

meta_conds = {
        "2016" : "$CMSSW_BASE/src/flashgg/MetaData/data/MetaConditions/Era2016_RR-17Jul2018_v1.json",
        "2017" : "$CMSSW_BASE/src/flashgg/MetaData/data/MetaConditions/Era2017_RR-31Mar2018_v1.json",
        "2018" : "$CMSSW_BASE/src/flashgg/MetaData/data/MetaConditions/Era2018_RR-17Sep2018_v1.json"
}

for coupling in couplings.keys():
    for year in years:
        couplings[coupling][year] = {}
        couplings[coupling][year]["outdir"] = "workspaces_" + coupling + "_" + year + "_" + args.tag + "_" + args.date + "/"
        os.system("mkdir -p %s" %  couplings[coupling][year]["outdir"] )
        os.system("cp workspaceStd.py %s" %  couplings[coupling][year]["outdir"] )

tarfile = "package_" + args.tag + ".tar.gz"

if os.path.isfile(tarfile) and not args.save_tar:
    os.system("rm %s" % tarfile)

# Make 1 tarball
cmssw_ver = [x for x in os.path.expandvars("$CMSSW_BASE").split("/") if "CMSSW" in x][0]
if not args.dry_run and not args.save_tar:
    print "Using CMSSW version: ", cmssw_ver
    os.system("XZ_OPT='-9e -T24' tar -Jc --exclude='.git' --exclude='my*.root' --exclude='*.tar*' --exclude='merged_ntuple*.root' --exclude='*.out' --exclude='*.err' --exclude='*.log' --exclude '*Taggers/data/DNN_models/*' --exclude '*Taggers/data/*dijet*' --exclude '*76X*.root' --exclude '*76X*.xml' --exclude '*Taggers/data/DNN_models/*' --exclude '*Taggers/data/HHTagger/*' --exclude '*Taggers/data/ttHKiller/*' -f %s ../../../../../%s/" % (tarfile, cmssw_ver))


productions = {}
procs = args.procs.split(",")
# Send tarball to multiple locations
for coupling in couplings.keys():
    for year in years:
        couplings[coupling][year]["stage_dir"] = "/hadoop/cms/store/user/smay/FCNC/" + couplings[coupling][year]["outdir"]
        os.system("mkdir -p %s" % couplings[coupling][year]["stage_dir"])
        if not args.dry_run and not args.save_tar:
            os.system("cp %s %s" % (tarfile, couplings[coupling][year]["stage_dir"] + "package.tar.gz"))
            os.system("hadoop fs -setrep -R 30 %s" % (couplings[coupling][year]["stage_dir"].replace("/hadoop","") + "package.tar.gz"))

        couplings[coupling][year]["datasets"] = {}
        for proc in procs:
            if proc == "fcnc":
                continue
            with open("wsJSONs/legacy_runII_v1_%s_%s.json" % (year, proc), "r") as f_in:
                datasets_to_add = json.load(f_in)
            for process, info in datasets_to_add["processes"].iteritems():
                if isinstance(info[0], list):
                    info = info[0]
                production = choose_production(info[0])                
                couplings[coupling][year]["datasets"][process] = { "sample" : info[0].split("/")[1], "catalog_name" : info[0], "cmdLine" : datasets_to_add["cmdLine"] , "production" : production}
                productions[info[0].split("/")[1]] = production
             

catalogs = ["Era2016_RR-17Jul2018_v2", "Era2017_RR-31Mar2018_v2", "Era2018_RR-17Sep2018_v2"]
samples = get_samples_from_catalogs(catalogs, productions)

datasets = {
        "2016" : {},
        "2017" : {},
        "2018" : {}
        }

for dataset in datasets.keys():
    for proc in procs:
        datasets[dataset][proc] = {}

for dataset in [dir for dir in glob.glob(args.fcnc_wildcard) if "RunIIFall17MiniAODv2" in dir]:
    sample = dataset.split("/")[-1]
    datasets["2017"]["fcnc"][sample] = { "input_loc" : dataset }

for year in years:
    for proc in procs:
        if proc == "fcnc":
            continue
        with open("wsJSONs/legacy_runII_v1_%s_%s.json" % (year, proc), "r") as f_in:
            datasets_to_add = json.load(f_in)
        for sample, info in datasets_to_add["processes"].iteritems():
            if isinstance(info[0], list):
                info = info[0]
            sample_name = info[0].split("/")[1]
            catalog_name = info[0]
            cmdLine = datasets_to_add["cmdLine"]
            if sample_name in samples.keys():
                xs = samples[sample_name]["xs"]
                scale1fb = samples[sample_name][year]["scale1fb"]
                production_names = []
                for production in samples[sample_name][year].keys():
                    if production == "scale1fb":
                        continue
                    else:
                        production_names.append(production)
                if len(production_names) > 1 and not ("DoubleEG" in sample or "EGamma" in sample):
                    print "More than one possible production to pick for sample %s" % sample_name
                if len(production_names) == 0:
                    print "Did not identify production for sample %s, year %s" % (sample_name, year)
                    continue
                files = []
                nEvents = 0
                weights = 0

                for production in production_names:
                    files += samples[sample_name][year][production]["files"]
                    nEvents += samples[sample_name][year][production]["nevents"]
                    weights += samples[sample_name][year][production]["weights"]

                datasets[year][proc][sample] = {
                        "sample_name" : sample_name,
                        "catalog_name" : catalog_name,
                        "cmdLine" : cmdLine,
                        "xs" : xs,
                        "scale1fb" : scale1fb,
                        "production_names" : production_names,
                        "files" : files,
                        "nEvents" : "nEvents",
                        "weights" : weights,
                        "args" : cmdLine + " dataset=%s" % catalog_name + " processId=%s" % sample,
                        "proc" : proc,
                        "year" : year,
                        }

            else:
                print "Could not find sample %s in flashgg catalog!!!" % sample_name


metadata = {
        "job_tag" : args.tag,
        "cmssw_ver" : cmssw_ver,
}

fpo = {
        "data" : 25,
        "sig" : 1,
        "fcnc" : 1,
        "zh" : 1,
}

command_list = []


for coupling, info in couplings.iteritems():
    if coupling == "Hut":
        coupling_selection = "fcncHutTagsOnly=True"
    elif coupling == "Hct":
        coupling_selection = "fcncHctTagsOnly=True"
    syst_selection = "doSystematics=%s" % syst
    for proc in procs:
        proc_name = "FCNC" if proc == "fcnc" else proc
        for year in years:
            if len(datasets[year].keys()) == 0:
                continue
            metadata_ = copy.deepcopy(metadata)
            metadata_["job_tag"] += "_" + coupling

            fcnc_executable = "%s/executable_%s.sh" % (info[year]["outdir"], proc)
            os.system("cp fcnc_executable.sh %s" % fcnc_executable)
            os.system("sed -i 's@MY_DIR@%s@g' %s" % (info[year]["outdir"], fcnc_executable))

            if proc == "fcnc":
                command = "$CMSSW_BASE/src/flashgg/Systematics/test/workspaceStd.py doHTXS=True useAAA=True doPdfWeights=False processId=%s %s %s %s %s %s" % (proc_name, coupling_selection, syst_selection, "outputFile=" + couplings[coupling][year]["outdir"] + "/output_FCNC_USER.root", "jobId=-1", "metaConditions=" + meta_conds[year]) 
            else:
                command = "$CMSSW_BASE/src/flashgg/Systematics/test/workspaceStd.py doHTXS=True doPdfWeights=False %s %s %s %s" % (coupling_selection, syst_selection, "outputFile=" + couplings[coupling][year]["outdir"] + "/output_FCNC_USER.root", "jobId=-1")
            os.system("sed -i 's@COMMAND@%s@g' %s" % (command, fcnc_executable))

            metadata_["executable"] = fcnc_executable
            metadata_["tarfile"] = info[year]["outdir"] + "package.tar.gz"
            metadata_["hadoop_path"] = info[year]["stage_dir"].replace("/hadoop/cms/store/user/smay/","")
            metadata_["args"] = ""
            metadata_["fpo"] = fpo[proc]

            metadata_file = info[year]["outdir"] + "%s_metadata.json" % proc
            with open(metadata_file, "w") as f_out:
                json.dump(metadata_, f_out, sort_keys=True, indent=4)

            datasets_file = info[year]["outdir"] + "%s_datasets.json" % proc
            with open(datasets_file, "w") as f_out:
                json.dump(datasets[year][proc], f_out, sort_keys=True, indent=4)

            function = "python submit_jobs.py --proc %s --metadata %s --datasets %s" % (proc, metadata_file, datasets_file)

            couplings[coupling][year]["command_" + proc] = function
            command_list.append(function)

# Bookkeeping
with open("ws_jobs_summary_%s.json" % args.tag, "w") as f_out:
    json.dump(couplings, f_out, sort_keys=True, indent=4)
with open("ws_datasets_%s.json" % args.tag, "w") as f_out:
    json.dump(datasets, f_out, sort_keys=True, indent=4)

# Submit all
if args.dry_run:
    for command in command_list:
        print "\n\nCOMMAND"
        print command
else:
    #parallel_utils.submit_jobs(command_list, 1, True, True)
    parallel_utils.submit_jobs(command_list, len(command_list), True, True)

import sys, os
import json
import datetime

d = datetime.datetime.today()

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--tag", help = "tag to identify this set of babies", type=str)
parser.add_argument("--years", help = "which year(s) to submit jobs for", type=str, default="2016,2017,2018")
parser.add_argument("--date", help = "date to distinguish this submission", type=str, default=d.strftime("%d-%B-%Y"))
parser.add_argument("--procs", help = "csv list of types of processes to submit", type=str, default = "sig,data,zh")
parser.add_argument("--couplings", help = "csv list of couplings", type=str, default = "Hct,Hut")
parser.add_argument("--dry_run", help = "don't actually make tarball or submit", action="store_true")
args = parser.parse_args()

import parallel_utils

couplings = {}
for coupling in args.couplings.split(","):
    couplings[coupling] = {}

years = args.years.split(",")

for coupling in couplings.keys():
    for year in years:
        couplings[coupling][year] = {}
        couplings[coupling][year]["outdir"] = "workspaces_" + coupling + "_" + year + "_" + args.tag + "_" + args.date + "/"
        os.system("mkdir -p %s" %  couplings[coupling][year]["outdir"] )
        os.system("cp workspaceStd.py %s" %  couplings[coupling][year]["outdir"] )

tarfile = "package_" + args.tag + ".tar.gz"

if os.path.isfile(tarfile):
    os.system("rm %s" % tarfile)

# Make 1 tarball
if not args.dry_run:
    cmssw_ver = [x for x in os.path.expandvars("$CMSSW_BASE").split("/") if "CMSSW" in x][0]
    print "Using CMSSW version: ", cmssw_ver
    os.system("XZ_OPT='-9e -T24' tar -Jc --exclude='.git' --exclude='my*.root' --exclude='*.tar*' --exclude='merged_ntuple*.root' --exclude='*.out' --exclude='*.err' --exclude='*.log' --exclude '*Taggers/data/DNN_models/*' --exclude '*Taggers/data/*dijet*' --exclude '*76X*.root' --exclude '*76X*.xml' --exclude '*Taggers/data/DNN_models/*' --exclude '*Taggers/data/HHTagger/*' --exclude '*Taggers/data/ttHKiller/*' -f %s ../../../../../%s/" % (tarfile, cmssw_ver))

# Send tarball to multiple locations
for coupling in couplings.keys():
    for year in years:
        couplings[coupling][year]["stage_dir"] = "/hadoop/cms/store/user/smay/FCNC/" + couplings[coupling][year]["outdir"]
        os.system("mkdir -p %s" % couplings[coupling][year]["stage_dir"])
        if not args.dry_run:
            os.system("cp %s %s" % (tarfile, couplings[coupling][year]["stage_dir"] + "package.tar.gz"))
            os.system("hadoop fs -setrep -R 30 %s" % (couplings[coupling][year]["stage_dir"].replace("/hadoop","") + "package.tar.gz"))

command_list = []
procs = args.procs.split(",")

for coupling, info in couplings.iteritems():
    for proc in procs:
        for year in years:
            if coupling == "Hut":
                coupling_selection = "fcncHutTagsOnly=True"
            elif coupling == "Hct":
                coupling_selection = "fcncHctTagsOnly=True"
            syst_selection = "" if proc == "data" else "doSystematics=True"
            command = "fggRunJobs.py --load wsJSONs/legacy_runII_v1_%s_%s.json -d %s workspaceStd.py -n 300 --no-copy-proxy -b htcondor --stage-to %s -q workday doHTXS=True %s %s" % (year, proc, info[year]["outdir"], "gsiftp://gftp.t2.ucsd.edu" + info[year]["stage_dir"], coupling_selection, syst_selection)
            couplings[coupling][year]["command_" + proc] = command
            command_list.append(command)

# Bookkeeping
with open("ws_jobs_summary_%s.json" % args.tag, "w") as f_out:
    json.dump(couplings, f_out, sort_keys=True, indent=4)

# Submit all
if args.dry_run:
    for command in command_list:
        print command
else:
    parallel_utils.submit_jobs(command_list, len(command_list), False)

import os

os.system("python $CMSSW_BASE/src/flashgg/Systematics/scripts/hadd_all.py")

current_dir = os.popen("pwd").read()
if "bbhth" in current_dir:
  os.system("source $CMSSW_BASE/src/flashgg/Systematics/scripts/hadd_bbhth.sh")

else:
  os.system("python $CMSSW_BASE/src/flashgg/Systematics/scripts/split_all.py")
  os.system("python $CMSSW_BASE/src/flashgg/Systematics/scripts/HTXS_mv_irrelevant.py")
  os.system("rm *_FWDH.root")
  os.system("python $CMSSW_BASE/src/flashgg/Systematics/scripts/rename.py")

import sys, os
import time
import json

from metis.Sample import DirectorySample
from metis.CondorTask import CondorTask
from metis.StatsParser import StatsParser

def submit_fcnc_workspaces(datasets, metadata):
    total_summary = {}
    while True:
        all_complete = False
        for dataset, info in datasets.items():
            sample = DirectorySample( dataset = dataset, location = info["input_loc"])
            task = CondorTask(
                    sample = sample,
                    open_dataset = False,
                    flush = True,
                    files_per_output = metadata["fpo"],
                    output_name = "output_FCNC_USER.root",
                    tag = metadata["job_tag"],
                    cmssw_version = metadata["cmssw_ver"],
                    executable = metadata["executable"],
                    tarfile = "dummy_tar",
                    condor_submit_params = { "sites" : "T2_US_UCSD" },
                    special_dir = metadata["hadoop_path"],
                    arguments = metadata["args"]
                    )
            task.process()
            if not task.complete():
                all_complete = False
            total_summary[dataset] = task.get_task_summary()
        #StatsParser(data = total_summary, webdir = "~/public_html/dump/FCNC_Workspaces_" + metadata["job_tag"] + "/").do()
        os.system("chmod -R 755 ~/public_html/dump/FCNC_Workspaces/")
        if all_complete:
            print("\n Job={} finished\n".format(metadata["job_tag"]))
            break
        print("Sleeping 300 seconds ...")
        time.sleep(300)

    return

if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("--metadata", type=json.loads)
    p.add_argument("--datasets", type=json.loads)
    arguments = p.parse_args()

    submit_fcnc_workspaces(arguments.datasets, arguments.metadata)

# NB this command is specific to the configuration at IC and is not gaurenteed elsewhere
fggRunJobs.py --load big_test_jobs_25ns.json -d big_test_jobs_1_0_0_WS0 -x cmsRun MicroAODtoWorkspace.py maxEvents=-1 -n 500 -q heplong.q -D -P useAAA=1 --no-use-tarball 
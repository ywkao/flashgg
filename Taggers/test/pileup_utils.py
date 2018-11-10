import sys, os
import glob
import json

import numpy as np

import FWCore.ParameterSet.Config as cms

def create_sample_map():
  datasets = {}
  input_json = "datasets_RunIIFall17.json"
  with open(input_json, "r") as f_in:
    datasets = json.load(f_in)

  samples_to_reweight = []
  for key, info in datasets["processes"].iteritems():
    if key == "data":
      continue
    for sample in info["datasets"]:
      samples_to_reweight.append(str(sample.split("/")[1]))

  pu_mix_files = glob.glob("MetaData/python/PU_MixFiles_2017_miniaodv2_310/*.py")

  sample_map = {}
  for sample in samples_to_reweight:
    sample_map[sample] = {}
    sample_map[sample]["n_match"] = 0
    for file in pu_mix_files:
      if "PSWeights" in file:
        continue
      if sample in file:
        sample_map[sample]["pu_mix"] = file
        sample_map[sample]["n_match"] += 1
  return sample_map

def match_file_to_sample(filename, sample_map):
  for sample, info in sample_map.iteritems():
    if sample in filename:
      print "Matched %s to %s" % (filename, sample)
      return sample
  print "Did not match file to sample"
  return ""


def read_file(file):
  with open(file, 'r') as f_in:
    return f_in.readlines()

def grep_info_from_file(file):
  lines = read_file(file)

  for i in range(len(lines)):
    if "probValue" in lines[i]:
      values = lines[i+1]
  values = values.split(",")
  values = [float(i) for i in values] 
 
  return cms.vdouble(values)      

def set_pileup_reweighting(obj, filename):
  sample_map = create_sample_map()
  sample = match_file_to_sample(filename, sample_map)
  if sample == "":
    obj.puReWeight = False
  else:
    mcPu = grep_info_from_file(sample_map[sample]["pu_mix"])
    print mcPu

    obj.puReWeight = True
    obj.puBins = cms.vdouble(map(float, cms.vint32(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99)))
    obj.mcPu = mcPu
    obj.useTruePu = cms.bool(True)

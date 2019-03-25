import sys, os

n_events = -1 
file_names = "DoubleEG"
#file_names = sys.argv[1].replace("/hadoop", "file:/hadoop").split(",")
#meta_conditions = sys.argv[2]

import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
import FWCore.ParameterSet.VarParsing as VarParsing
from flashgg.Systematics.SystematicDumperDefaultVariables import minimalVariables,minimalHistograms,minimalNonSignalVariables,systematicVariables
from flashgg.Systematics.SystematicDumperDefaultVariables import minimalVariablesHTXS,systematicVariablesHTXS
import os
from flashgg.MetaData.MetaConditionsReader import *

# SYSTEMATICS SECTION
process = cms.Process("FLASHggTag")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32( 1000 )


#process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring("file:/hadoop/cms/store/user/smay/ttH/MicroAOD/RunII/DoubleEG_Run2016B-17Jul2018_ver2-v1_MINIAOD_RunII/test_skim_1.root"))
process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring("file:/hadoop/cms/store/user/smay/ttH/MicroAOD/RunII/DoubleEG_Run2017C-31Mar2018-v1_MINIAOD_RunII/test_skim_2.root"))

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(n_events))


systlabels = [""]
phosystlabels = []
metsystlabels = []
jetsystlabels = []
elesystlabels = []
musystlabels = []

from flashgg.MetaData.MetaConditionsReader import *

# set default options if needed
#metaConditions = MetaConditionsReader("MetaData/data/MetaConditions/Era2016_RR-17Jul2018_v1.json")
metaConditions = MetaConditionsReader("MetaData/data/MetaConditions/Era2017_RR-31Mar2018_v1.json")
#metaConditions = MetaConditionsReader("MetaData/data/MetaConditions/Era2018_RR-17Sep2018_v1.json")

ISDATA = False
if "DoubleEG" in file_names[0] or "EGamma" in file_names[0]:
  ISDATA = True

ISSIGNAL = False
if "ttH" in file_names[0]:
  ISSIGNAL = True


### Global Tag
from Configuration.AlCa.GlobalTag import GlobalTag
if ISDATA:
    process.GlobalTag.globaltag = str(metaConditions['globalTags']['data'])
else:
    process.GlobalTag.globaltag = str(metaConditions['globalTags']['MC'])

from flashgg.Systematics.SystematicsCustomize import *
class customizer:
  def __init__(self, **kwargs):
    self.metaConditions = kwargs.get("metaConditions")
    self.processType = kwargs.get("processType")

customize = customizer(metaConditions = metaConditions, processType = "data" if ISDATA else "mc")
jetSystematicsInputTags = createStandardSystematicsProducers(process , customize)
modifyTagSequenceForSystematics(process,jetSystematicsInputTags)

# needed for 0th vertex from microAOD
process.load("flashgg/MicroAOD/flashggDiPhotons_cfi")
process.flashggDiPhotons.whichVertex = cms.uint32(0)
process.flashggDiPhotons.useZerothVertexFromMicro = cms.bool(True)

# Remove unneeded tags
process.flashggTagSequence.remove(process.flashggTTHDiLeptonTag)
process.flashggTagSequence.remove(process.flashggVBFTag)
process.flashggTagSequence.remove(process.flashggVHMetTag)
process.flashggTagSequence.remove(process.flashggWHLeptonicTag)
process.flashggTagSequence.remove(process.flashggZHLeptonicTag)
process.flashggTagSequence.remove(process.flashggVHLeptonicLooseTag)
process.flashggTagSequence.remove(process.flashggVHHadronicTag)
process.flashggTagSequence.remove(process.flashggUntagged)
process.flashggTagSequence.remove(process.flashggVBFMVA)
process.flashggTagSequence.remove(process.flashggVBFDiPhoDiJetMVA)

process.flashggTagSorter.TagPriorityRanges = cms.VPSet(   
        cms.PSet(TagName = cms.InputTag('flashggTTHLeptonicTag')),
        cms.PSet(TagName = cms.InputTag('flashggTTHHadronicTag')) )

newvpset = cms.VPSet()
for pset in process.flashggDiPhotonSystematics.SystMethods:
    if not pset.Label.value().count("FracRVNvtxWeight") :
	print  pset.Label.value()
	newvpset += [pset]
#from flashgg.Systematics.flashggDiPhotonSystematics_cfi import PixelSeedWeight
#newvpset += [ PixelSeedWeight ]

process.flashggDiPhotonSystematics.SystMethods = newvpset


# Or use the official tool instead
useEGMTools(process)

if ISSIGNAL:
  print "Signal MC, so adding systematics and dZ"
  customizeSystematicsForSignal(process)
elif ISDATA:
  print "Data, so turn off all shifts and systematics, with some exceptions"
  variablesToUse = minimalNonSignalVariables
  customizeSystematicsForData(process)
else:
  print "Background MC, so store mgg and central only"
  variablesToUse = minimalNonSignalVariables
  customizeSystematicsForBackground(process)

printSystematicInfo(process)

from HLTrigger.HLTfilters.hltHighLevel_cfi import hltHighLevel
hlt_paths = []
for dset in metaConditions["TriggerPaths"]:
  hlt_paths.extend(metaConditions["TriggerPaths"][dset])
for i in range(len(hlt_paths)):
  hlt_paths[i] = hlt_paths[i].encode("ascii")
process.hltHighLevel= hltHighLevel.clone(HLTPaths = cms.vstring(hlt_paths))

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

# ee bad supercluster filter on data
process.load('RecoMET.METFilters.eeBadScFilter_cfi')
process.eeBadScFilter.EERecHitSource = cms.InputTag("reducedEgamma","reducedEERecHits") # Saved MicroAOD Collection (data only)
process.dataRequirements = cms.Sequence()
if ISDATA: 
        process.dataRequirements += process.hltHighLevel
        process.dataRequirements += process.eeBadScFilter


process.genFilter = cms.Sequence()

process.load("flashgg/Taggers/flashggTagSequence_cfi")
process.load("flashgg/Taggers/flashggTagTester_cfi")

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("merged_ntuple.root"))
## TAGS DUMPERS ##
from flashgg.Taggers.tagsDumpers_cfi import *

process.tthLeptonicTagDumper = createTagDumper("TTHLeptonicTag")
process.tthHadronicTagDumper = createTagDumper("TTHHadronicTag")

process.tthLeptonicTagDumper.dumpTrees = True
process.tthHadronicTagDumper.dumpTrees = True

## define categories and associated objects to dump
import flashgg.Taggers.dumperConfigTools as cfgTools

dipho_variables=["dipho_sumpt      := diPhoton.sumPt",
                 "dipho_cosphi     := abs(cos(diPhoton.leadingPhoton.phi - diPhoton.subLeadingPhoton.phi))",
                 "mass             := diPhoton.mass",
                 "leadPt           := diPhoton.leadingPhoton.pt",
                 "leadEt           := diPhoton.leadingPhoton.et",
		 "leadEnergy       := diPhoton.leadingPhoton.energy",
                 "leadEta          := diPhoton.leadingPhoton.eta",
                 "leadPhi          := diPhoton.leadingPhoton.phi",
                 "lead_sieie       := diPhoton.leadingPhoton.sigmaIetaIeta",
                 "lead_hoe         := diPhoton.leadingPhoton.hadronicOverEm",
                 "lead_sigmaEoE    := diPhoton.leadingPhoton.sigEOverE",
                 "lead_ptoM        := diPhoton.leadingPhoton.pt/diPhoton.mass",
                 "leadR9           := diPhoton.leadingPhoton.full5x5_r9"]

## TTH leptonic ##
cfgTools.addCategories(process.tthLeptonicTagDumper,
                       ## categories definition
                       [("all","1",0)
                    ],
                       ## variables to be dumped in trees/datasets. Same variables for all categories
                       variables=dipho_variables+
                       ["n_ele    := electrons.size"
		],
                       ## histograms
                       histograms=[]
)

## TTH hadronic ##
cfgTools.addCategories(process.tthHadronicTagDumper,
                       ## categories definition
                       [("all","1",0)
                    ],
                       ## variables to be dumped in trees/datasets. Same variables for all categories
                       variables=dipho_variables+
                       ["n_bjets  := nBMedium",
                        "n_jets   := jetVector.size"
		],
                       ## histograms
                       histograms=[]
)

process.p = cms.Path(process.dataRequirements*
                         #process.genFilter*
                         process.flashggDiPhotons* # needed for 0th vertex from microAOD
                         process.flashggUpdatedIdMVADiPhotons*
                         #process.flashggDiPhotonSystematics*
                         process.flashggMetSystematics*
                         process.flashggMuonSystematics*process.flashggElectronSystematics*
                         (process.flashggUnpackedJets*process.jetSystematicsSequence)*
                         (process.flashggTagSequence*process.systematicsTagSequences)*
			 process.flashggTagSequence*
                         process.flashggTagTester*
		         (process.tthLeptonicTagDumper
                          +process.tthHadronicTagDumper))

# Comment
import sys

import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
from flashgg.Systematics.SystematicDumperDefaultVariables import minimalVariables,minimalHistograms,minimalNonSignalVariables,systematicVariables

process = cms.Process("FLASHggTag")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#process.GlobalTag.globaltag = 'POSTLS170_V5::All'
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc')
nEvents = -1
if len(sys.argv) > 3:
  nEvents = int(sys.argv[3])

print('Running over %d events' % nEvents)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(nEvents) )
#process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32( 100 )

MUON_ID = "Tight" #["Tight", "Medium" , "Loose", "Soft", "HighPt", "MediumPrompt", "TrkHighPt"]
MUON_ISO = "LooseRel" #{ LooseID : ["LooseRel"],MediumID:["LooseRel", "TightRel"] , TrkHighPtID:["LooseRelTk", "TightRelTk"], TightIDandIPCut:["LooseRel", "TightRel"], HighPtIDandIPCut:["LooseRelTk", "TightRelTk"] }
from flashgg.Systematics.SystematicsCustomize import *
jetSystematicsInputTags = createStandardSystematicsProducers(process)
modifyTagSequenceForSystematics(process,jetSystematicsInputTags)

# Uncomment the following if you notice you have a memory leak
# This is a lightweight tool to digg further
#process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck",
#                                        ignoreTotal = cms.untracked.int32(1),
#                                        monitorPssAndPrivate = cms.untracked.bool(True)
#                                       )

#process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring("file:/hadoop/cms/store/user/smay/tth/ttHJetToGG_M100_13TeV_amcatnloFXFX_madspin_pythia8/RunIISummer16-2_4_1-25ns_Moriond17-2_4_1-v0-RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/170119_110959/0000/myMicroAODOutputFile_1.root"))
fileNames = sys.argv[2]
fileNames = fileNames.replace('/hadoop','file:/hadoop')
fileNames = fileNames.split(",")
process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring(fileNames))
#process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring("file:/hadoop/cms/store/user/bemarsh/flashgg/MicroAOD_skim/2016_skim_v2/DoubleEG/ReMiniAOD-03Feb2017-2_5_4-2_5_1-v0-Run2016B-03Feb2017_ver2-v2/microAOD_1.root", "file:/hadoop/cms/store/user/bemarsh/flashgg/MicroAOD_skim/2016_skim_v2/DoubleEG/ReMiniAOD-03Feb2017-2_5_4-2_5_1-v0-Run2016B-03Feb2017_ver2-v2/microAOD_2.root"))

# Add triggers/filters from workspaceStd.py L444-463
from HLTrigger.HLTfilters.hltHighLevel_cfi import hltHighLevel
process.hltHighLevel= hltHighLevel.clone(HLTPaths = cms.vstring("HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90_v*",
#                                                                "HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_DoublePixelVeto_Mass55_v1",
#                                                                "HLT_Diphoton30EB_18EB_R9Id_OR_IsoCaloId_AND_HE_R9Id_DoublePixelVeto_Mass55_v1"
                                                                ))

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

# ee bad supercluster filter on data
process.load('RecoMET.METFilters.eeBadScFilter_cfi')
process.eeBadScFilter.EERecHitSource = cms.InputTag("reducedEgamma","reducedEERecHits") # Saved MicroAOD Collection (data only)
# Bad Muon filter LOADS WRONG IN 8_0_28, FIX LATER
#process.load('RecoMET.METFilters.badGlobalMuonTaggersMiniAOD_cff')
#process.badGlobalMuonTaggerMAOD.muons = cms.InputTag("flashggSelectedMuons")
#process.cloneGlobalMuonTaggerMAOD.muons = cms.InputTag("flashggSelectedMuons")
process.dataRequirements = cms.Sequence()

#print('Process ID: %s' % customize.processId)
#if customize.processId == "Data":
if "DoubleEG" in fileNames[0]:
        print('Adding HLT and filter requirements')
        process.dataRequirements += process.hltHighLevel
        process.dataRequirements += process.eeBadScFilter

process.load("flashgg/Taggers/flashggTagSequence_cfi")
process.load("flashgg/Taggers/flashggTagTester_cfi")

print process.flashggTagSequence

useEGMTools(process)

if "ttH" in fileNames[0]:
  print "Signal MC, so adding systematics and dZ"
  customizeSystematicsForSignal(process)
elif "DoubleEG" in fileNames[0]:
  print "Data, so turn off all shifts and systematics, with some exceptions"
  variablesToUse = minimalNonSignalVariables
  customizeSystematicsForData(process)
else:
  print "Background MC, so store mgg and central only"
  variablesToUse = minimalNonSignalVariables
  customizeSystematicsForBackground(process)

printSystematicInfo(process)

#if "ttH" in fileNames[0]:
#  customizeSystematicsForSignal(process)
#elif "DoubleEG" in fileNames[0]:
#  customizeSystematicsForData(process)
#else:
#  customizeSystematicsForBackground(process):

# For debugging
switchToUnPreselected = False
switchToFinal = False
switchToPuppi = False
switchToReadOld = False
assert(not switchToUnPreselected or not switchToFinal)
assert(not switchToReadOld or not switchToUnPreselected)
assert(not switchToReadOld or not switchToFinal)

if switchToReadOld:
    from PhysicsTools.PatAlgos.tools.helpers import massSearchReplaceAnyInputTag
    massSearchReplaceAnyInputTag(process.flashggTagSequence,cms.InputTag("flashggPreselectedDiPhotons"),cms.InputTag("flashggDiPhotonsWithAddedDz"))
    process.flashggDiPhotonsWithAddedDz = cms.EDProducer('FlashggDiPhotonGenZProducer',
                                                 DiPhotonTag=cms.InputTag('flashggPreselectedDiPhotons'),
                                                 GenParticleTag=cms.InputTag( "flashggPrunedGenParticles" ))
    process.flashggNewPreselectedDiPhotons = cms.Sequence(process.flashggPreselectedDiPhotons*process.flashggDiPhotonsWithAddedDz)
    process.flashggTagSequence.replace(process.flashggPreselectedDiPhotons,process.flashggNewPreselectedDiPhotons)
    process.source.fileNames=cms.untracked.vstring("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/cmshgg/sethzenz/flashgg/RunIISpring15-ReMiniAOD-BetaV7-25ns/Spring15BetaV7/GluGluHToGG_M-125_13TeV_powheg_pythia8/RunIISpring15-ReMiniAOD-BetaV7-25ns-Spring15BetaV7-v0-RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/151021_152108/0000/myMicroAODOutputFile_2.root")
    print process.flashggTagSequence

if switchToUnPreselected:
    from PhysicsTools.PatAlgos.tools.helpers import massSearchReplaceAnyInputTag
    massSearchReplaceAnyInputTag(process.flashggTagSequence,cms.InputTag("flashggPreselectedDiPhotons"),cms.InputTag("flashggDiPhotons"))

if switchToFinal:
    from flashgg.MicroAOD.flashggFinalEGamma_cfi import flashggFinalEGamma
    from PhysicsTools.PatAlgos.tools.helpers import massSearchReplaceAnyInputTag
    massSearchReplaceAnyInputTag(process.flashggTagSequence,cms.InputTag("flashggPreselectedDiPhotons"),cms.InputTag("flashggFinalEGamma",flashggFinalEGamma.DiPhotonCollectionName.value()))

if switchToPuppi:
    process.flashggUnpackedJets.JetsTag = cms.InputTag("flashggFinalPuppiJets")


from flashgg.Taggers.flashggTagOutputCommands_cff import tagDefaultOutputCommand

#process.out = cms.OutputModule("PoolOutputModule", fileName = cms.untracked.string('myTagOutputFile.root'),
#                               outputCommands = tagDefaultOutputCommand			       
#                               )

process.p = cms.Path(process.dataRequirements*
		     process.flashggUpdatedIdMVADiPhotons*
		     process.flashggDiPhotonSystematics*
		     process.flashggMetSystematics*
		     process.flashggMuonSystematics*process.flashggElectronSystematics*
		     (process.flashggUnpackedJets*process.jetSystematicsSequence)*
		     (process.flashggTagSequence*process.systematicsTagSequences)*
		     process.flashggTagSequence*
		     process.flashggTagTester)

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
                 "leadEta          := diPhoton.leadingPhoton.eta",
                 "leadPhi          := diPhoton.leadingPhoton.phi",
                 "lead_sieie       := diPhoton.leadingPhoton.sigmaIetaIeta",
                 "lead_hoe         := diPhoton.leadingPhoton.hadronicOverEm",
                 "lead_sigmaEoE    := diPhoton.leadingPhoton.sigEOverE",
                 "lead_ptoM        := diPhoton.leadingPhoton.pt/diPhoton.mass",
                 "leadR9           := diPhoton.leadingPhoton.full5x5_r9",
                 "leadGenMatch     := diPhoton.leadingPhoton.genMatchType",
                 "subleadPt        := diPhoton.subLeadingPhoton.pt",
                 "subleadEt        := diPhoton.subLeadingPhoton.et",
                 "subleadEta       := diPhoton.subLeadingPhoton.eta",
                 "subleadPhi       := diPhoton.subLeadingPhoton.phi",
                 "sublead_sieie    := diPhoton.subLeadingPhoton.sigmaIetaIeta",
                 "sublead_hoe      := diPhoton.subLeadingPhoton.hadronicOverEm",
                 "sublead_sigmaEoE := diPhoton.subLeadingPhoton.sigEOverE",
                 "sublead_ptoM     := diPhoton.subLeadingPhoton.pt/diPhoton.mass",
                 "subleadR9        := diPhoton.subLeadingPhoton.full5x5_r9",
                 "subleadGenMatch  := diPhoton.subLeadingPhoton.genMatchType",
                 "leadIDMVA        := diPhoton.leadingView.phoIdMvaWrtChosenVtx",
                 "subleadIDMVA     := diPhoton.subLeadingView.phoIdMvaWrtChosenVtx",
                 "dipho_rapidity   := diPhoton.rapidity",
                 "vertex_idx       := diPhoton.vertexIndex",
                 "nGoodEls         := nGoodEls",
                 "nGoodElsFromTau  := nGoodElsFromTau",
                 "nGoodMus         := nGoodMus",
                 "nGoodMusFromTau  := nGoodMusFromTau",
                 "nGoodTaus        := nGoodTaus",
                 "diphoMVARes      := diphoMVARes",
                 "leadPassEVeto    := diPhoton.leadingPhoton.passElectronVeto",
                 "subleadPassEVeto := diPhoton.subLeadingPhoton.passElectronVeto",
                 "leadPixelSeed    := diPhoton.leadingPhoton.hasPixelSeed",
                 "subleadPixelSeed := diPhoton.subLeadingPhoton.hasPixelSeed",
                 "nb_loose         := nBLoose",
                 "nb_medium        := nBMedium",
                 "nb_tight         := nBTight",
                 "lead_photon_type := leadPhotonType",
                 "sublead_photon_type := subleadPhotonType"
]

## TTH leptonic ##
cfgTools.addCategories(process.tthLeptonicTagDumper,
                       ## categories definition
                       [("all","1",0)
                    ],
                       ## variables to be dumped in trees/datasets. Same variables for all categories
                       variables=dipho_variables+
                       ["n_ele    := electrons.size",
                        "ele1_pt  := ?(electrons.size>0)? electrons.at(0).pt : -1",
                        "ele2_pt  := ?(electrons.size>1)? electrons.at(1).pt : -1",
                        "n_muons  := muons.size",
                        "muon1_pt := ?(muons.size>0)? muons.at(0).pt : -1",
                        "muon2_pt := ?(muons.size>1)? muons.at(1).pt : -1",
                        "n_bjets  := bJets.size",
                        "n_jets   := jets.size",
                        "bjet1_pt := bJets.at(0).pt",
                        "bjet2_pt := ?bJets.size>1? bJets.at(1).pt : -1",
                        "bjet1_csv:= bJets.at(0).bDiscriminator('pfCombinedInclusiveSecondaryVertexV2BJetTags')",
                        "bjet2_csv:= ?bJets.size>1? bJets.at(1).bDiscriminator('pfCombinedInclusiveSecondaryVertexV2BJetTags') : -1",
                        "Mjj      := ?(jets.size>1)?"
                        +"sqrt((jets.at(0).energy+jets.at(1).energy)^2-(jets.at(0).px+jets.at(1).px)^2-(jets.at(0).py+jets.at(1).py)^2-(jets.at(0).pz+jets.at(1).pz)^2)"
                        +": -1",
                        "MetPt  := MetPt",
                        "MetPhi := MetPhi",
                        "mT     := MT",
                        "jet_pt1                :=  ? jets.size()>0 ? jets[0].pt() : -100 ",
                        "jet_eta1               :=  ? jets.size()>0 ? jets[0].eta() : -100 ",
                        "jet_phi1               :=  ? jets.size()>0 ? jets[0].phi() : -100 ",
                        "jet_bdiscriminant1     :=  ? jets.size()>0 ? jets[0].bDiscriminator('pfCombinedInclusiveSecondaryVertexV2BJetTags') : -100",
                        "jet_pt2                :=  ? jets.size()>1 ? jets[1].pt() : -100 ",
                        "jet_eta2               :=  ? jets.size()>1 ? jets[1].eta() : -100 ",
                        "jet_phi2               :=  ? jets.size()>1 ? jets[1].phi() : -100 ",
                        "jet_bdiscriminant2     :=  ? jets.size()>1 ? jets[1].bDiscriminator('pfCombinedInclusiveSecondaryVertexV2BJetTags') : -100",
                        "jet_pt3                :=  ? jets.size()>2 ? jets[2].pt() : -100 ",
                        "jet_eta3               :=  ? jets.size()>2 ? jets[2].eta() : -100 ",
                        "jet_phi3               :=  ? jets.size()>2 ? jets[2].phi() : -100 ",
                        "jet_bdiscriminant3     :=  ? jets.size()>2 ? jets[2].bDiscriminator('pfCombinedInclusiveSecondaryVertexV2BJetTags') : -100 ",
                        "jet_pt4                :=  ? jets.size()>3 ? jets[3].pt() : -100 ",
                        "jet_eta4               :=  ? jets.size()>3 ? jets[3].eta() : -100 ",
                        "jet_phi4               :=  ? jets.size()>3 ? jets[3].phi() : -100 ",
                        "jet_bdiscriminant4     :=  ? jets.size()>3 ? jets[3].bDiscriminator('pfCombinedInclusiveSecondaryVertexV2BJetTags') : -100 ",
                        "jet_pt5                :=  ? jets.size()>4 ? jets[4].pt() : -100 ",
                        "jet_eta5               :=  ? jets.size()>4 ? jets[4].eta() : -100 ",
                        "jet_phi5               :=  ? jets.size()>4 ? jets[4].phi() : -100 ",
                        "jet_bdiscriminant5     :=  ? jets.size()>4 ? jets[4].bDiscriminator('pfCombinedInclusiveSecondaryVertexV2BJetTags') : -100 ",
                        "jet_pt6                :=  ? jets.size()>5 ? jets[5].pt() : -100 ",
                        "jet_eta6               :=  ? jets.size()>5 ? jets[5].eta() : -100 ",
                        "jet_phi6               :=  ? jets.size()>5 ? jets[5].phi() : -100 ",
                        "jet_bdiscriminant6     :=  ? jets.size()>5 ? jets[5].bDiscriminator('pfCombinedInclusiveSecondaryVertexV2BJetTags') : -100 ",
                        "jet_pt7                :=  ? jets.size()>6 ? jets[6].pt() : -100 ",
                        "jet_eta7               :=  ? jets.size()>6 ? jets[6].eta() : -100 ",
                        "jet_phi7               :=  ? jets.size()>6 ? jets[6].phi() : -100 ",
                        "jet_bdiscriminant7     :=  ? jets.size()>6 ? jets[6].bDiscriminator('pfCombinedInclusiveSecondaryVertexV2BJetTags') : -100 ",
                        "jet_pt8                :=  ? jets.size()>7 ? jets[7].pt() : -100 ",
                        "jet_eta8               :=  ? jets.size()>7 ? jets[7].eta() : -100 ",
                        "jet_phi8               :=  ? jets.size()>7 ? jets[7].phi() : -100 ",
                        "jet_bdiscriminant8     :=  ? jets.size()>7 ? jets[7].bDiscriminator('pfCombinedInclusiveSecondaryVertexV2BJetTags') : -100 ",
"jet_pt9                :=  ? jets.size()>8 ? jets[8].pt() : -100 ",
                        "jet_eta9               :=  ? jets.size()>8 ? jets[8].eta() : -100 ",
                        "jet_phi9               :=  ? jets.size()>8 ? jets[8].phi() : -100 ",
                        "jet_bdiscriminant9     :=  ? jets.size()>8 ? jets[8].bDiscriminator('pfCombinedInclusiveSecondaryVertexV2BJetTags') : -100 ",
                        "jet_pt10                :=  ? jets.size()>9 ? jets[9].pt() : -100 ",
                        "jet_eta10               :=  ? jets.size()>9 ? jets[9].eta() : -100 ",
                        "jet_phi10               :=  ? jets.size()>9 ? jets[9].phi() : -100 ",
                        "jet_bdiscriminant10     :=  ? jets.size()>9 ? jets[9].bDiscriminator('pfCombinedInclusiveSecondaryVertexV2BJetTags') : -100 ",
                        "jet_pt11                :=  ? jets.size()>10 ? jets[10].pt() : -100 ",
                        "jet_eta11               :=  ? jets.size()>10 ? jets[10].eta() : -100 ",
                        "jet_phi11               :=  ? jets.size()>10 ? jets[10].phi() : -100 ",
                        "jet_bdiscriminant11     :=  ? jets.size()>10 ? jets[10].bDiscriminator('pfCombinedInclusiveSecondaryVertexV2BJetTags') : -100 ",
                        "jet_pt12                :=  ? jets.size()>11 ? jets[11].pt() : -100 ",
                        "jet_eta12               :=  ? jets.size()>11 ? jets[11].eta() : -100 ",
                        "jet_phi12               :=  ? jets.size()>11 ? jets[11].phi() : -100 ",
                        "jet_bdiscriminant12     :=  ? jets.size()>11 ? jets[11].bDiscriminator('pfCombinedInclusiveSecondaryVertexV2BJetTags') : -100 ",
                        "jet_pt13                :=  ? jets.size()>12 ? jets[12].pt() : -100 ",
                        "jet_eta13               :=  ? jets.size()>12 ? jets[12].eta() : -100 ",
                        "jet_phi13               :=  ? jets.size()>12 ? jets[12].phi() : -100 ",
                        "jet_bdiscriminant13     :=  ? jets.size()>12 ? jets[12].bDiscriminator('pfCombinedInclusiveSecondaryVertexV2BJetTags') : -100 ",
                        "jet_pt14                :=  ? jets.size()>13 ? jets[13].pt() : -100 ",
                        "jet_eta14               :=  ? jets.size()>13 ? jets[13].eta() : -100 ",
                        "jet_phi14               :=  ? jets.size()>13 ? jets[13].phi() : -100 ",
                        "jet_bdiscriminant14     :=  ? jets.size()>13 ? jets[13].bDiscriminator('pfCombinedInclusiveSecondaryVertexV2BJetTags') : -100 ",
                        "jet_pt15                :=  ? jets.size()>14 ? jets[14].pt() : -100 ",
                        "jet_eta15               :=  ? jets.size()>14 ? jets[14].eta() : -100 ",
                        "jet_phi15               :=  ? jets.size()>14 ? jets[14].phi() : -100 ",
                        "jet_bdiscriminant15     :=  ? jets.size()>14 ? jets[14].bDiscriminator('pfCombinedInclusiveSecondaryVertexV2BJetTags') : -100 ",


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
                        "n_jets   := jetVector.size",
                        "bjet1_pt := bJetVector.at(0).pt",
                        "bjet2_pt := ?nBMedium>1? bJetVector.at(1).pt : -1",
                        "MetPt  := MetPt",
                        "MetPhi := MetPhi",
                        "jet1_pt  := ?(jetVector.size>0)? jetVector.at(0).pt : -1",
                        "jet2_pt  := ?(jetVector.size>1)? jetVector.at(1).pt : -1",
                        "jet3_pt  := ?(jetVector.size>2)? jetVector.at(2).pt : -1",
                        "jet4_pt  := ?(jetVector.size>3)? jetVector.at(3).pt : -1",
                        "jet5_pt  := ?(jetVector.size>4)? jetVector.at(4).pt : -1",
                        "jet6_pt  := ?(jetVector.size>5)? jetVector.at(5).pt : -1",
                        "jet7_pt  := ?(jetVector.size>6)? jetVector.at(6).pt : -1",
                        "jet8_pt  := ?(jetVector.size>7)? jetVector.at(7).pt : -1",
                        "jet9_pt  := ?(jetVector.size>8)? jetVector.at(8).pt : -1",
                        "jet10_pt  := ?(jetVector.size>9)? jetVector.at(9).pt : -1",
                        "jet11_pt  := ?(jetVector.size>10)? jetVector.at(10).pt : -1",
                        "jet12_pt  := ?(jetVector.size>11)? jetVector.at(11).pt : -1",
                        "jet13_pt  := ?(jetVector.size>12)? jetVector.at(12).pt : -1",
                        "jet14_pt  := ?(jetVector.size>13)? jetVector.at(13).pt : -1",
                        "jet15_pt  := ?(jetVector.size>14)? jetVector.at(14).pt : -1",
                        "jet1_eta  := ?(jetVector.size>0)? jetVector.at(0).eta : -1",
                        "jet2_eta  := ?(jetVector.size>1)? jetVector.at(1).eta : -1",
                        "jet3_eta  := ?(jetVector.size>2)? jetVector.at(2).eta : -1",
                        "jet4_eta  := ?(jetVector.size>3)? jetVector.at(3).eta : -1",
                        "jet5_eta  := ?(jetVector.size>4)? jetVector.at(4).eta : -1",
                        "jet6_eta  := ?(jetVector.size>5)? jetVector.at(5).eta : -1",
                        "jet7_eta  := ?(jetVector.size>6)? jetVector.at(6).eta : -1",
                        "jet8_eta  := ?(jetVector.size>7)? jetVector.at(7).eta : -1",
                        "jet9_eta  := ?(jetVector.size>8)? jetVector.at(8).eta : -1",
                        "jet10_eta  := ?(jetVector.size>9)? jetVector.at(9).eta : -1",
                        "jet11_eta  := ?(jetVector.size>10)? jetVector.at(10).eta : -1",
                        "jet12_eta  := ?(jetVector.size>11)? jetVector.at(11).eta : -1",
                        "jet13_eta  := ?(jetVector.size>12)? jetVector.at(12).eta : -1",
                        "jet14_eta  := ?(jetVector.size>13)? jetVector.at(13).eta : -1",
                        "jet15_eta  := ?(jetVector.size>14)? jetVector.at(14).eta : -1",
                        "jet1_phi  := ?(jetVector.size>0)? jetVector.at(0).phi : -1",
                        "jet2_phi  := ?(jetVector.size>1)? jetVector.at(1).phi : -1",
                        "jet3_phi  := ?(jetVector.size>2)? jetVector.at(2).phi : -1",
                        "jet4_phi  := ?(jetVector.size>3)? jetVector.at(3).phi : -1",
                        "jet5_phi  := ?(jetVector.size>4)? jetVector.at(4).phi : -1",
                        "jet6_phi  := ?(jetVector.size>5)? jetVector.at(5).phi : -1",
                        "jet7_phi  := ?(jetVector.size>6)? jetVector.at(6).phi : -1",
                        "jet8_phi  := ?(jetVector.size>7)? jetVector.at(7).phi : -1",
                        "jet9_phi  := ?(jetVector.size>8)? jetVector.at(8).phi : -1",
                        "jet10_phi  := ?(jetVector.size>9)? jetVector.at(9).phi : -1",
                        "jet11_phi  := ?(jetVector.size>10)? jetVector.at(10).phi : -1",
                        "jet12_phi  := ?(jetVector.size>11)? jetVector.at(11).phi : -1",
                        "jet13_phi  := ?(jetVector.size>12)? jetVector.at(12).phi : -1",
                        "jet14_phi  := ?(jetVector.size>13)? jetVector.at(13).phi : -1",
                        "jet15_phi  := ?(jetVector.size>14)? jetVector.at(14).phi : -1",
                        "bjet1_csv:= maxBTagVal",
                        "bjet2_csv:= secondMaxBTagVal",
                        "tthMVA   := tthMvaRes"
                    ],
                       ## histograms
                       histograms=[]
)

process.p1 = cms.Path(process.dataRequirements*
		      (process.tthLeptonicTagDumper
                      +process.tthHadronicTagDumper))
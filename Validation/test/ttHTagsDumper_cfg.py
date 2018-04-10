#!/usr/bin/env cmsRun

import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
from FWCore.ParameterSet.VarParsing import VarParsing

from flashgg.MetaData.samples_utils import SamplesManager

## CMD LINE OPTIONS ##
options = VarParsing('analysis')

# maxEvents is the max number of events processed of each file, not globally
options.maxEvents = -1
options.inputFiles = "file:myTagOutputFile.root"
options.outputFile = "merged_ntuple.root"
options.parseArguments()

## I/O SETUP ##
process = cms.Process("ValidationTagsDumper")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32( 1 )
process.source = cms.Source ("PoolSource",
                             fileNames = cms.untracked.vstring(options.inputFiles))#,

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string(options.outputFile))

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
                 "subleadPt        := diPhoton.subLeadingPhoton.pt",
                 "subleadEt        := diPhoton.subLeadingPhoton.et",
                 "subleadEta       := diPhoton.subLeadingPhoton.eta",
                 "subleadPhi       := diPhoton.subLeadingPhoton.phi",
                 "sublead_sieie    := diPhoton.subLeadingPhoton.sigmaIetaIeta",
                 "sublead_hoe      := diPhoton.subLeadingPhoton.hadronicOverEm",
                 "sublead_sigmaEoE := diPhoton.subLeadingPhoton.sigEOverE",
                 "sublead_ptoM     := diPhoton.subLeadingPhoton.pt/diPhoton.mass",
                 "subleadR9        := diPhoton.subLeadingPhoton.full5x5_r9",
                 "leadIDMVA        := diPhoton.leadingView.phoIdMvaWrtChosenVtx",
                 "subleadIDMVA     := diPhoton.subLeadingView.phoIdMvaWrtChosenVtx",
		 "dipho_rapidity   := diPhoton.rapidity",
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
			"bjet1_csv:= maxBTagVal",
			"bjet2_csv:= secondMaxBTagVal",
			"tthMVA   := tthMvaRes"
                    ],
                       ## histograms
                       histograms=[]
)

process.p1 = cms.Path(process.tthLeptonicTagDumper
                      +process.tthHadronicTagDumper)

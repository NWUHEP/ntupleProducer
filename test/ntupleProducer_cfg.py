import os
import FWCore.ParameterSet.Config as cms
from RecoEgamma.PhotonIdentification.isolationCalculator_cfi import *

process = cms.Process("PAT")

# real data or MC?
isRealData = False

# global tag
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

if (isRealData):
    process.GlobalTag.globaltag = 'GR_R_44_V14::All' 
else:
    process.GlobalTag.globaltag = 'START44_V13::All'

# tau reconstruction configuration
process.load("RecoTauTag.Configuration.RecoPFTauTag_cff")

# jet energy corrections
process.load('JetMETCorrections.Configuration.DefaultJEC_cff')

# pat sequences
process.load("PhysicsTools.PatAlgos.patSequences_cff")

from PhysicsTools.PatAlgos.patEventContent_cff import patEventContent
process.out = cms.OutputModule("PoolOutputModule",
                               fileName = cms.untracked.string('/tmp/patTuple.root'),
                               SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('p') ),
                               outputCommands = cms.untracked.vstring('drop *', *patEventContent )
                               )

from Higgs.ntupleProducer.PatSequences_cff import addPatSequence
addPatSequence(process, isRealData, addPhotons = True)

##------------------------------------------------------------------
## see PhysicsTools/PatExamples/test/analyzePatTau_fromAOD_cfg.py
## switch the tau algorithm (AA))
#from PhysicsTools.PatAlgos.tools.tauTools import *
#switchToPFTauHPS(process)

# global options
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'ERROR'
process.MessageLogger.cerr.FwkReport.reportEvery = 100000000

'''
process.MessageLogger.categories = cms.untracked.vstring('FwkJob', 'FwkReport', 'FwkSummary', 'Root_NoDictionary', 'DataNotAvailable', 'HLTConfigData')
process.MessageLogger.destinations = cms.untracked.vstring('myOutput')
process.MessageLogger.myOutput = cms.untracked.PSet(
                FwkJob              = cms.untracked.PSet(limit = cms.untracked.int32(0)),
                FwkReport           = cms.untracked.PSet(limit = cms.untracked.int32(0)),
                FwkSummary          = cms.untracked.PSet(limit = cms.untracked.int32(0)),
                Root_NoDictionary   = cms.untracked.PSet(limit = cms.untracked.int32(0)),                    
                DataNotAvailable    = cms.untracked.PSet(limit = cms.untracked.int32(0)),                    
                HLTConfigData       = cms.untracked.PSet(limit = cms.untracked.int32(0))
                )

process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(False),
                                     SkipEvent = cms.untracked.vstring('ProductNotFound')
                                    )
'''

# event source  
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(10))
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
       '/store/data/Fall11/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/AODSIM/PU_S6-START44_V5-v1/0000/0030ADBC-C409-E111-B1E7-E0CB4E553666.root',
         #'/store/data/Run2011B/DoubleMu/AOD/16Jan2012-v1/0000/A0914A57-C344-E111-8687-001A928116D0.root'
         #'/store/user/stoyan/MC/H135toZG_500k/RECO_v1/H135toZG_7TeV_START44_V5_RAW2DIGI_RECO_PU_file9J_1_1_xse.root'
)
)

#process.load('gravitonData')


# event counters
process.startCounter = cms.EDProducer("EventCountProducer")
process.endCounter = process.startCounter.clone()

#process.ntuplePath = cms.Path(process.startCounter * process.patDefaultSequence * process.endCounter)

# configure output
from Higgs.ntupleProducer.OutputConfiguration_cff import configureOutput
configureOutput(process)

print '\n\nNow run the ntuplizer...\n\n'

### TFile service!
process.TFileService = cms.Service('TFileService',
                                   fileName = cms.string('nuTuple.root')
                                   )

### ntuple producer
process.ntupleProducer   = cms.EDAnalyzer('ntupleProducer',

  photonIsoCalcTag  =    cms.PSet(isolationSumsCalculator),

  JetTag            =    cms.untracked.InputTag('selectedPatJetsPFlow'),
  JPTTag            =    cms.untracked.InputTag('JetPlusTrackZSPCorJetAntiKt5'),
  GenJetTag         =    cms.untracked.InputTag('ak5GenJets'),
  METTag            =    cms.untracked.InputTag('patMETsPFlow'),
  METNoPUTag        =    cms.untracked.InputTag('patMETsPFlowNoPileup'),
  ElectronTag       =    cms.untracked.InputTag('selectedPatElectronsPFlow'),
  MuonTag           =    cms.untracked.InputTag('selectedPatMuonsPFlow'),
  #MuonTag           =    cms.untracked.InputTag('muons'),
  PhotonTag         =    cms.untracked.InputTag('patPhotons'),
  TauTag            =    cms.untracked.InputTag('selectedPatTausPFlow'),
  PrimaryVtxTag     =    cms.untracked.InputTag('offlinePrimaryVertices'),
  rhoCorrTag        =    cms.untracked.InputTag('kt6PFJetsPFlow', 'rho', 'PAT'),

  saveJets          =    cms.untracked.bool(True),
  saveElectrons     =    cms.untracked.bool(True),
  saveMuons         =    cms.untracked.bool(True),
  saveTaus          =    cms.untracked.bool(True),
  savePhotons       =    cms.untracked.bool(True),
  saveMET           =    cms.untracked.bool(True),
  saveGenJets       =    cms.untracked.bool(True),

  ecalFilterTag     =    cms.untracked.InputTag("BE1214","anomalousECALVariables"),
  hcalFilterTag     =    cms.untracked.InputTag("HBHENoiseFilterResultProducer","HBHENoiseFilterResult"),

  hltName           =    cms.untracked.string("HLT"),
  triggers          =    cms.untracked.vstring(
                                               "HLT_Mu8_v",
                                               "HLT_Mu15_v",
                                               "HLT_Mu8_Jet40_v",
                                               "HLT_Mu13_Mu8_v",
                                               "HLT_Mu17_Mu8_v",
                                               "HLT_DoubleMu3_v",
                                               "HLT_DoubleMu6_v",
                                               "HLT_DoubleMu7_v",
                                               "",
                                               "",
                                               "",
                                               "",
                                               "HLT_Ele8_CaloIdL_CaloIsoVL_v",
                                               "HLT_Ele17_CaloIdL_CaloIsoVL_v",
                                               "HLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v",
                                               "HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v",
                                               "HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v",
                                               "HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v",
                                               "",
                                               "",
                                               "",
                                               "",
                                               "HLT_Photon20_CaloIdVL_IsoL_v",
                                               "HLT_Photon20_CaloIdVL_v",
                                               "HLT_Photon30_CaloIdVL_IsoL_v",
                                               "HLT_Photon30_CaloIdVL_v",
                                               "HLT_Photon50_CaloIdVL_IsoL_v",
                                               "HLT_Photon50_CaloIdVL_v",
                                               "HLT_Photon75_CaloIdVL_IsoL_v",
                                               "HLT_Photon75_CaloIdVL_v",
                                               "HLT_Photon90_CaloIdVL_IsoL_v",
                                               "HLT_Photon90_CaloIdVL_v",
                                               "HLT_Photon135_v",
                                               "",
                                               "",
                                               "",
                                               "",
                                               "HLT_Mu17_Ele8_CaloIdL_v",
                                               "HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v",
                                               "HLT_Mu3_Ele8_CaloIdT_CaloIsoVL_v",
                                               "HLT_Mu8_Ele17_CaloIdL_v",
                                               "HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v",
                                               "",
                                               "",
                                               "",
                                               "",
                                               "HLT_Ele18_CaloIdVT_TrkIdT_MediumIsoPFTau20_v",
                                               "HLT_Ele20_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_MediumIsoPFTau20_v",
                                               "HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_MediumIsoPFTau25_v",
                                               "HLT_HT400_DoubleIsoPFTau10_Trk3_PFMHT50_v",
                                               "HLT_IsoMu15_eta2p1_LooseIsoPFTau20_v", 
                                               "HLT_IsoMu15_eta2p1_MediumIsoPFTau20_v", 
                                               "HLT_IsoMu15_eta2p1_TightIsoPFTau20_v", 
                                               "HLT_Mu15_LooseIsoPFTau15_v"
                                               "",
                                               "",
                                               "",
                                               ""
)
)

process.ntuplePath = cms.Path(process.PFTau
        * process.patDefaultSequence 
        * process.ntupleProducer
        )

#process.outpath = cms.EndPath(process.out)

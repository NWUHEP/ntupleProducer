from PhysicsTools.PatAlgos.patTemplate_cfg import *

process.load("RecoTauTag.Configuration.RecoPFTauTag_cff")

# real data or MC?
isRealData = False

# global tag
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

if (isRealData):
    process.GlobalTag.globaltag = 'GR_R_42_V20::All' 
else:
    process.GlobalTag.globaltag = 'START44_V12::All'

# jet energy corrections
process.load('JetMETCorrections.Configuration.DefaultJEC_cff')

# global options
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(False), SkipEvent = cms.untracked.vstring('ProductNotFound'))
process.MessageLogger.cerr.FwkReport.reportEvery = 5000

# event source  
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(100))
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
         #'/store/mc/Fall11/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/AODSIM/PU_S6_START42_V14B-v1/0001/A6A53A52-4EF6-E011-8524-90E6BA0D09AD.root'
         '/store/user/stoyan/MC/H135toZG_500k/RECO_v1/H135toZG_7TeV_START44_V5_RAW2DIGI_RECO_PU_file9J_1_1_xse.root'
)
)

# pat sequences
from Higgs.ntupleProducer.PatSequences_cff import addPatSequence
addPatSequence(process, isRealData, addPhotons = True)

##------------------------------------------------------------------
## see PhysicsTools/PatExamples/test/analyzePatTau_fromAOD_cfg.py
## switch the tau algorithm (AA))
#from PhysicsTools.PatAlgos.tools.tauTools import *
#switchToPFTauHPS(process)


# event counters
process.startCounter = cms.EDProducer("EventCountProducer")
process.endCounter = process.startCounter.clone()

#process.ntuplePath = cms.Path(process.startCounter * process.patDefaultSequence * process.endCounter)

# configure output
from Higgs.ntupleProducer.OutputConfiguration_cff import configureOutput
configureOutput(process)
#process.out.fileName = cms.untracked.string('test.root')

print '\n\nNow run the ntuplizer...\n\n'

### TFile service!
process.TFileService = cms.Service('TFileService',
                                   fileName = cms.string('nuTuple.root')
                                   )

### ntuple producer
process.ntupleProducer   = cms.EDAnalyzer('ntupleProducer',

  JetTag            =    cms.untracked.InputTag('selectedPatJetsPFlow'),
  GenJetTag         =    cms.untracked.InputTag('ak5GenJets'),
  METTag            =    cms.untracked.InputTag('patMETsPFlow'),
  ElectronTag       =    cms.untracked.InputTag('selectedPatElectronsPFlow'),
  MuonTag           =    cms.untracked.InputTag('selectedPatMuonsPFlow'),
  PhotonTag         =    cms.untracked.InputTag('patPhotons'),
#####  TauTag            =    cms.untracked.InputTag('shrinkingConePFTauProducer'),
# only the one below is in PAT -> check selections
  TauTag            =    cms.untracked.InputTag('selectedPatTausPFlow'),
  PrimaryVtxTag     =    cms.untracked.InputTag('offlinePrimaryVertices'),
  rhoCorrTag        =    cms.untracked.InputTag('kt6PFJetsIso', 'rho', 'ntuples'),

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
  #triggers          =    testTriggers
  triggers          =    cms.untracked.vstring(
                                               "HLT_Mu8_v",
                                               "HLT_Mu15_v",
                                               "HLT_Mu8_Jet40_v",
                                               "HLT_Mu13_Mu8_v",
                                               "HLT_Mu17_Mu8_v",
                                               "HLT_DoubleMu3_v",
                                               "HLT_DoubleMu6_v",
                                               "HLT_DoubleMu7_v",

                                               "HLT_Ele8_CaloIdL_CaloIsoVL_v",
                                               "HLT_Ele17_CaloIdL_CaloIsoVL_v",
                                               "HLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v",
                                               "HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v",
                                               "HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v",
                                               "HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v",

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

                                               "HLT_Mu17_Ele8_CaloIdL_v",
                                               "HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v",
                                               "HLT_Mu3_Ele8_CaloIdT_CaloIsoVL_v",
                                               "HLT_Mu8_Ele17_CaloIdL_v",
                                               "HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v"
)
)

#process.ntuplePath = cms.Path(process.PFTau * process.PFTauprocess.patDefaultSequence * process.ntupleProducer)
process.ntuplePath = cms.Path(process.PFTau * process.patDefaultSequence * process.ntupleProducer)

import FWCore.ParameterSet.Config as cms

process = cms.Process("ntuples") 

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'ERROR'

process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(10000000)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(500))

process.load("Configuration.StandardSequences.Services_cff")
process.load('Configuration.StandardSequences.GeometryExtended_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

### Conditions tags
process.GlobalTag.globaltag = 'GR_R_42_V20::All' 
#process.GlobalTag.globaltag = 'START42_V13::All' 


### HCAL noise filter
# Hcal Noise filter
process.load('CommonTools/RecoAlgos/HBHENoiseFilterResultProducer_cfi')
process.HBHENoiseFilterResultProducer.minNumIsolatedNoiseChannels = cms.int32(9999)
process.HBHENoiseFilterResultProducer.minIsolatedNoiseSumE = cms.double(9999)
process.HBHENoiseFilterResultProducer.minIsolatedNoiseSumEt = cms.double(9999)

### Ecal noise filter
process.load('PhysicsTools/EcalAnomalousEventFilter/ecalanomalouseventfilter_cfi')
process.EcalAnomalousEventFilter.FilterAlgo= cms.untracked.string("TuningMode")
process.EcalAnomalousEventFilter.cutBoundEnergyDeadCellsEB=cms.untracked.double(10)
process.EcalAnomalousEventFilter.cutBoundEnergyDeadCellsEE=cms.untracked.double(10)
process.EcalAnomalousEventFilter.cutBoundEnergyGapEB=cms.untracked.double(100)
process.EcalAnomalousEventFilter.cutBoundEnergyGapEE=cms.untracked.double(100)
process.EcalAnomalousEventFilter.enableGap=cms.untracked.bool(False)

process.BE1214 = process.EcalAnomalousEventFilter.clone()
process.BE1214.limitDeadCellToChannelStatusEB = cms.vint32(12,14)
process.BE1214.limitDeadCellToChannelStatusEE = cms.vint32(12,14)
process.load('JetMETAnalysis.ecalDeadCellTools.RA2TPfilter_cff')

ecalDead = cms.Sequence(process.BE1214)

### Jet correction services
process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
process.load("CondCore.DBCommon.CondDBCommon_cfi")

### Extra jet collection for L1FastJet corrections
process.load("RecoJets.Configuration.RecoPFJets_cff")
process.kt6PFJets.doRhoFastjet = True
process.kt6PFJets.Rho_EtaMax   = cms.double(4.4)
process.kt6PFJets.rParam       = cms.double(0.6)

process.kt6PFJetsIso = process.kt6PFJets.clone()
process.kt6PFJetsIso.doRhoFastjet = True
process.kt6PFJetsIso.Rho_EtaMax   = cms.double(2.5) 
process.ak5PFJets.doAreaFastjet   = True
process.ak5PFJets.Rho_EtaMax      = cms.double(4.4)

### To get b-tags from ak5PFJets
process.load('RecoJets.JetAssociationProducers.ak5JTA_cff')
process.ak5JetTracksAssociatorAtVertex.jets = cms.InputTag("ak5PFJetsL1FastL2L3")
process.ak5JetTracksAssociatorAtCaloFace.jets = cms.InputTag("ak5PFJetsL1FastL2L3")
process.ak5JetExtender.jets = cms.InputTag("ak5PFJetsL1FastL2L3")

### GenJet flavor matching
process.load("PhysicsTools.JetMCAlgos.CaloJetsMCFlavour_cfi")

# Flavour byReference 
process.GenJetbyRef = cms.EDProducer("JetPartonMatcher", 
                                     jets = cms.InputTag("ak5GenJets"), 
                                     coneSizeToAssociate = cms.double(0.3), 
                                     partons = cms.InputTag("myPartons") 
                                     ) 
# Flavour byValue PhysDef 
process.GenJetbyValPhys = cms.EDProducer("JetFlavourIdentifier", 
                                         srcByReference = cms.InputTag("GenJetbyRef"), 
                                         physicsDefinition = cms.bool(True), 
                                         leptonInfo = cms.bool(True) 
                                         ) 
# Flavour byValue AlgoDef 
process.GenJetbyValAlgo = cms.EDProducer("JetFlavourIdentifier", 
                                         srcByReference = cms.InputTag("GenJetbyRef"), 
                                         physicsDefinition = cms.bool(False), 
                                         leptonInfo = cms.bool(True) 
                                         ) 
process.GenJetFlavour = cms.Sequence(process.GenJetbyRef*process.GenJetbyValPhys*process.GenJetbyValAlgo)


process.JetbyRef = cms.EDProducer("JetPartonMatcher",
                                  jets = cms.InputTag("ak5PFJets"),
                                  coneSizeToAssociate = cms.double(0.3),
                                  partons = cms.InputTag("myPartons")
                                  )
# Flavour byValue PhysDef
process.JetbyValPhys = cms.EDProducer("JetFlavourIdentifier",
                                      srcByReference = cms.InputTag("JetbyRef"),
                                      physicsDefinition = cms.bool(True),
                                      leptonInfo = cms.bool(True)
                                      )
# Flavour byValue AlgoDef
process.JetbyValAlgo = cms.EDProducer("JetFlavourIdentifier",
                                      srcByReference = cms.InputTag("JetbyRef"),
                                      physicsDefinition = cms.bool(False),
                                      leptonInfo = cms.bool(True)
                                      )
process.JetFlavour = cms.Sequence(process.JetbyRef*process.JetbyValPhys*process.JetbyValAlgo)

## eleID map:
process.load("ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi")
process.simpleEleId60relIso = process.simpleCutBasedElectronID.clone()
process.simpleEleId60relIso.electronQuality = "60relIso"
process.simpleEleId70relIso = process.simpleCutBasedElectronID.clone()
process.simpleEleId70relIso.electronQuality = "70relIso"
process.simpleEleId80relIso = process.simpleCutBasedElectronID.clone()
process.simpleEleId80relIso.electronQuality = "80relIso"
process.simpleEleId85relIso = process.simpleCutBasedElectronID.clone()
process.simpleEleId85relIso.electronQuality = "85relIso"
process.simpleEleId90relIso = process.simpleCutBasedElectronID.clone()
process.simpleEleId90relIso.electronQuality = "90relIso"
process.simpleEleId95relIso = process.simpleCutBasedElectronID.clone()
process.simpleEleId95relIso.electronQuality = "95relIso"

### MET corrections
process.load("JetMETCorrections.Type1MET.pfMETCorrections_cff")
#process.pfJetMETcorr.jetCorrLabel = cms.string("ak5PFL1FastL2L3Residual")
process.metAnalysisSequence=cms.Sequence(process.producePFMETCorrections)

### Specify default triggers
#from UserCode.ntupleProducer.triggerTest_cfi import *

### Filter by primary vertices
process.goodVertices = cms.EDFilter("VertexSelector",
  src = cms.InputTag('offlinePrimaryVertices'),
  cut = cms.string('!isFake && isValid && ndof >= 4.0 &&	position.Rho < 2.0 && abs(z) < 24'),
  filter = cms.bool(True)
)

### Input files
#process.load('H145toZG')

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
         '/store/data/Run2011A/DoubleMu/RECO/PromptReco-v4/000/165/121/1A873C93-A381-E011-902F-0030487CD710.root'
#         '/store/mc/Summer11/DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia/AODSIM/PU_S4_START42_V11-v1/0001/AC9DFDE2-7CA8-E011-B14E-002590200948.root'
)
)

### TFile service!
process.TFileService = cms.Service('TFileService',
                                   fileName = cms.string('nuTuple.root')
                                   )

### ntuple producer
process.ntupleProducer   = cms.EDAnalyzer('ntupleProducer',

  JetTag            =    cms.untracked.InputTag('ak5PFJets'),
  GenJetTag         =    cms.untracked.InputTag('ak5GenJets'),
  METTag            =    cms.untracked.InputTag('pfMet'),
  ElectronTag       =    cms.untracked.InputTag('gsfElectrons'),
  MuonTag           =    cms.untracked.InputTag('muons'),
  PhotonTag         =    cms.untracked.InputTag('photons'),
  TauTag            =    cms.untracked.InputTag('shrinkingConePFTauProducer'),
  PrimaryVtxTag     =    cms.untracked.InputTag('offlinePrimaryVertices'),
  rhoCorrTag        =    cms.untracked.InputTag('kt6PFJetsIso', 'rho', 'ntuples'),

  saveJets          =    cms.untracked.bool(True),
  saveElectrons     =    cms.untracked.bool(True),
  saveMuons         =    cms.untracked.bool(True),
  saveTaus          =    cms.untracked.bool(True),
  savePhotons       =    cms.untracked.bool(True),
  saveMET           =    cms.untracked.bool(True),
  saveGenJets       =    cms.untracked.bool(False),

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

                                               "HLT_Mu17_Ele8_CaloIdL_v",
                                               "HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v",
                                               "HLT_Mu3_Ele8_CaloIdT_CaloIsoVL_v",
                                               "HLT_Mu8_Ele17_CaloIdL_v",
                                               "HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v"
)
)

### Let it run
cmsSeq = cms.Sequence(
        process.PFTau                    
      #* process.myPartons #<-- For genJet flavors, only in MC
      #* process.GenJetFlavour
      #* process.JetFlavour
      * process.simpleEleId60relIso
      * process.simpleEleId70relIso
      * process.simpleEleId80relIso
      * process.simpleEleId85relIso
      * process.simpleEleId90relIso
      * process.simpleEleId95relIso
      * process.kt6PFJets
      * process.kt6PFJetsIso
      * process.ak5PFJets
      * process.ak5PFJetsL1FastL2L3
      * process.metAnalysisSequence  
      * process.ak5JetTracksAssociatorAtVertex 
      * process.btagging
      * process.HBHENoiseFilterResultProducer
		)

process.p = cms.Path(cmsSeq * process.ntupleProducer)

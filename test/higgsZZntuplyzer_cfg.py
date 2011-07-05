import FWCore.ParameterSet.Config as cms

process = cms.Process("higgs") 

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("FWCore.MessageLogger.MessageLogger_cfi")
#process.MessageLogger.cerr.threshold = 'ERROR'

process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(100)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(500))

process.load("Configuration.StandardSequences.Services_cff")
process.load('Configuration.StandardSequences.GeometryExtended_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

### Conditions tags
process.GlobalTag.globaltag = 'GR_R_42_V14::All' 

### Input files
#process.load("Jet")

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
#        'file:/tmp/naodell/62CCADC4-EC7B-E011-AC6F-0015178C4994.root'
#       '/store/data/Run2010A/JetMET/RECO/Sep17ReReco_v2/0026/FE9C387C-24C8-DF11-8CAA-003048679274.root'
'/store/data/Run2011A/DoubleMu/RECO/PromptReco-v4/000/165/121/1A873C93-A381-E011-902F-0030487CD710.root'
)
)

### HCAL noise filter
# Hcal Noise filter
#process.load('CommonTools/RecoAlgos/HBHENoiseFilter_cfi')
process.load('CommonTools/RecoAlgos/HBHENoiseFilterResultProducer_cfi')
process.HBHENoiseFilterResultProducer.minNumIsolatedNoiseChannels = cms.int32(99999)
process.HBHENoiseFilterResultProducer.minIsolatedNoiseSumE = cms.double(99999)
process.HBHENoiseFilterResultProducer.minIsolatedNoiseSumEt = cms.double(99999)

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

### Select primary vertices
process.load("RecoVertex.PrimaryVertexProducer.OfflinePrimaryVerticesDA_cfi")
process.offlinePrimaryVerticesDA = process.offlinePrimaryVerticesDA.clone()
process.offlinePrimaryVerticesDA.useBeamConstraint = cms.bool(True)
process.offlinePrimaryVerticesDA.TkClusParameters.TkDAClusParameters.Tmin= cms.double(4.)
process.offlinePrimaryVerticesDA.TkClusParameters.TkDAClusParameters.vertexSize= cms.double(0.01)

### Jet correction services
process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
process.load("CondCore.DBCommon.CondDBCommon_cfi")

### Producing the pho values for isolation (as from Jacub's email)
process.load("RecoJets.Configuration.RecoPFJets_cff")
process.kt6PFJets.doRhoFastjet = True
process.kt6PFJets.Rho_EtaMax = cms.double(4.4)
process.kt6PFJets.rParam = cms.double(0.6) 

process.kt6PFJetsIso = process.kt6PFJets.clone()
process.kt6PFJetsIso.doRhoFastjet = True
process.kt6PFJetsIso.Rho_EtaMax = cms.double(2.5) # this is used for rho calculation

process.ak5PFJets = process.ak5PFJets.clone()
process.ak5PFJets.doAreaFastjet = True
process.ak5PFJets.Rho_EtaMax = cms.double(4.4)

### To get b-tags from ak5PFJets
process.load('RecoJets.JetAssociationProducers.ak5JTA_cff')
process.ak5JetTracksAssociatorAtVertex.jets = cms.InputTag("ak5PFJetsL1FastL2L3")
process.ak5JetTracksAssociatorAtCaloFace.jets = cms.InputTag("ak5PFJetsL1FastL2L3")
process.ak5JetExtender.jets = cms.InputTag("ak5PFJetsL1FastL2L3")

### eleID map:
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
from JetMETCorrections.Type1MET.MetType1Corrections_cff import metJESCorAK5PFJet
process.metJESCorPFAK5 = metJESCorAK5PFJet.clone()
process.metJESCorPFAK5.inputUncorJetsLabel = "ak5PFJets"
process.metJESCorPFAK5.metType = "PFMET"
process.metJESCorPFAK5.inputUncorMetLabel = "pfMet"
process.metJESCorPFAK5.useTypeII = False
process.metJESCorPFAK5.jetPTthreshold = cms.double(10.0)
#process.metJESCorPFAK5.corrector = cms.string('ak5PFL1FastL2L3')
process.metJESCorPFAK5.corrector = cms.string('ak5PFL2L3')

### ntuple producer
process.ntupleProducer   = cms.EDAnalyzer('ntupleProducer',
#  rootfilename      =    cms.untracked.string("nuTuple_Photon_Run2011A.root"),

  JetTag            =    cms.untracked.InputTag("ak5PFJets"),
  GenJetTag         =    cms.untracked.InputTag("ak5GenJets"),
  METTag            =    cms.untracked.InputTag("pfMet"),
  ElectronTag       =    cms.untracked.InputTag("gsfElectrons"),
  MuonTag           =    cms.untracked.InputTag("muons"),
  PhotonTag         =    cms.untracked.InputTag("photons"),
  TauTag            =    cms.untracked.InputTag("shrinkingConePFTauProducer"),
  PrimaryVtxTag     =    cms.untracked.InputTag("offlinePrimaryVerticesDA"),
  rhoCorrTag	     =    cms.untracked.InputTag("kt6PFJetsIso", "rho", "higgs"),
  saveJets          =    cms.untracked.bool(True),
  saveElectrons     =    cms.untracked.bool(True),
  saveMuons         =    cms.untracked.bool(True),
  saveTaus          =    cms.untracked.bool(False),
  savePhotons       =    cms.untracked.bool(True),
  saveMET           =    cms.untracked.bool(True),
  saveGenJets       =    cms.untracked.bool(True),
                                          
  ecalAnomalousFilterTag = cms.untracked.InputTag("BE1214","anomalousECALVariables"),
  hcalFilterTag          = cms.untracked.InputTag("HBHENoiseFilterResultProducer","HBHENoiseFilterResult"),

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

                                               "HLT_Ele8_CaloIdL_CaloIsoVL_v",
                                               "HLT_Ele17_CaloIdL_CaloIsoVL_v",
                                               "HLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v",
                                               "HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v",
                                               "HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v",

                                               "HLT_Photon20_CaloIdL_CaloIsoVL_v",
                                               "HLT_Photon30_CaloIdL_CaloIsoVL_v",
                                               "HLT_Photon50_CaloIdL_CaloIsoVL_v",
                                               "HLT_Photon75_CaloIdL_CaloIsoVL_v",
                                               "HLT_Photon90_CaloIdL_CaloIsoVL_v"
)
)

process.load("NWU/ntupleProducer/hzzSkim_cff")

process.demo0 = cms.EDAnalyzer('Dummy',)
process.demo1 = process.demo0.clone()
process.demo2 = process.demo0.clone()
process.demo3 = process.demo0.clone()

process.options = cms.untracked.PSet(  wantSummary = cms.untracked.bool(True)    )

process.TimerService = cms.Service("TimerService", useCPUtime = cms.untracked.bool(True))
process.pts = cms.EDFilter("PathTimerInserter")
process.PathTimerService = cms.Service("PathTimerService")

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string('out_higgs.root')
                                   )


### Let it run
cmsSeq = cms.Sequence(
    process.offlinePrimaryVerticesDA 
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
    * process.ak5JetTracksAssociatorAtVertex 
    * process.metJESCorPFAK5  # corrected Met type I
    * process.btagging
    * process.HBHENoiseFilterResultProducer
    * ecalDead  #comment this out if don't want ecal filter
    )

process.diSequence = cms.Path(process.demo0
                              *(process.goodHzzMuons + process.goodHzzElectrons)
                              *(process.diHzzMuons + process.diHzzElectrons + process.crossHzzLeptons ) )

process.diMuonFilter     = cms.Path(process.diHzzMuonsFilter      *process.demo1* cmsSeq *process.ntupleProducer)
process.diElectronFilter = cms.Path(process.diHzzElectronsFilter  *process.demo2* cmsSeq *process.ntupleProducer)
process.EleMuFilter      = cms.Path(process.crossHzzLeptonsFilter *process.demo3* cmsSeq *process.ntupleProducer)


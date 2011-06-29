import FWCore.ParameterSet.Config as cms

process = cms.Process("ntuples") 

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'ERROR'

process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(100)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100))

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
        'file:/tmp/naodell/FEB6499E-8C7B-E011-98DA-0018F3D09702.root'
#       '/store/data/Run2010A/JetMET/RECO/Sep17ReReco_v2/0026/FE9C387C-24C8-DF11-8CAA-003048679274.root'
)
)

### HCAL noise filter
# Hcal Noise filter
#process.load('CommonTools/RecoAlgos/HBHENoiseFilter_cfi')
process.load('CommonTools/RecoAlgos/HBHENoiseFilterResultProducer_cfi')
process.HBHENoiseFilterResultProducer.minRatio = cms.double(-999)
process.HBHENoiseFilterResultProducer.maxRatio = cms.double(999)
process.HBHENoiseFilterResultProducer.minHPDHits = cms.int32(17)
process.HBHENoiseFilterResultProducer.minRBXHits = cms.int32(999)
process.HBHENoiseFilterResultProducer.minHPDNoOtherHits = cms.int32(10)
process.HBHENoiseFilterResultProducer.minZeros = cms.int32(10)
process.HBHENoiseFilterResultProducer.minHighEHitTime = cms.double(-9999.0)
process.HBHENoiseFilterResultProducer.maxHighEHitTime = cms.double(9999.0)
process.HBHENoiseFilterResultProducer.maxRBXEMF = cms.double(-999.0)
process.HBHENoiseFilterResultProducer.minNumIsolatedNoiseChannels = cms.int32(9999)
process.HBHENoiseFilterResultProducer.minIsolatedNoiseSumE = cms.double(9999)
process.HBHENoiseFilterResultProducer.minIsolatedNoiseSumEt = cms.double(9999)
process.HBHENoiseFilterResultProducer.useTS4TS5 = cms.bool(True)

### Ecal noise filter
# process.load('PhysicsTools/EcalAnomalousEventFilter/ecalanomalouseventfilter_cfi')
# process.EcalAnomalousEventFilter.FilterAlgo= cms.untracked.string("TuningMode")
# process.EcalAnomalousEventFilter.cutBoundEnergyDeadCellsEB=cms.untracked.double(10)
# process.EcalAnomalousEventFilter.cutBoundEnergyDeadCellsEE=cms.untracked.double(10)
# process.EcalAnomalousEventFilter.cutBoundEnergyGapEB=cms.untracked.double(100)
# process.EcalAnomalousEventFilter.cutBoundEnergyGapEE=cms.untracked.double(100)
# process.EcalAnomalousEventFilter.enableGap=cms.untracked.bool(False)
# 
# process.BE1214 = process.EcalAnomalousEventFilter.clone()
# process.BE1214.limitDeadCellToChannelStatusEB = cms.vint32(12,14)
# process.BE1214.limitDeadCellToChannelStatusEE = cms.vint32(12,14)
# process.load('JetMETAnalysis.ecalDeadCellTools.RA2TPfilter_cff')
# 
# ecalDead = cms.Sequence(process.BE1214)

### Select primary vertices
process.load("RecoVertex.PrimaryVertexProducer.OfflinePrimaryVerticesDA_cfi")
process.offlinePrimaryVerticesDAWithBS = process.offlinePrimaryVerticesDA.clone()
process.offlinePrimaryVerticesDAWithBS.useBeamConstraint = cms.bool(True)
process.offlinePrimaryVerticesDAWithBS.TkClusParameters.TkDAClusParameters.Tmin= cms.double(4.)
process.offlinePrimaryVerticesDAWithBS.TkClusParameters.TkDAClusParameters.vertexSize= cms.double(0.01)

### Jet correction services
process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
process.load("CondCore.DBCommon.CondDBCommon_cfi")

### Extra jet collection for L1FastJet corrections
process.load("RecoJets.Configuration.RecoPFJets_cff")
process.kt6PFJets.doRhoFastjet = True
process.kt6PFJets.Rho_EtaMax = cms.double(2.5)
process.kt6PFJets.rParam = cms.double(0.6)

process.ak5PFJets.doAreaFastjet = True
process.ak5PFJets.Rho_EtaMax = cms.double(2.5)

### To get b-tags from ak5PFJets
process.load('RecoJets.JetAssociationProducers.ak5JTA_cff')
process.ak5JetTracksAssociatorAtVertex.jets = cms.InputTag("ak5PFJetsL1FastL2L3")
process.ak5JetTracksAssociatorAtCaloFace.jets = cms.InputTag("ak5PFJetsL1FastL2L3")
process.ak5JetExtender.jets = cms.InputTag("ak5PFJetsL1FastL2L3")

### eleID map:
process.load("ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi")
process.simpleEleId80relIso = process.simpleCutBasedElectronID.clone()
process.simpleEleId80relIso.electronQuality = "80relIso"

### MET corrections
from JetMETCorrections.Type1MET.MetType1Corrections_cff import metJESCorAK5PFJet
process.metJESCorPFAK5 = metJESCorAK5PFJet.clone()
process.metJESCorPFAK5.inputUncorJetsLabel = "ak5PFJets"
process.metJESCorPFAK5.metType = "PFMET"
process.metJESCorPFAK5.inputUncorMetLabel = "pfMet"
process.metJESCorPFAK5.useTypeII = False
process.metJESCorPFAK5.jetPTthreshold = cms.double(10.0)
process.metJESCorPFAK5.corrector = cms.string('ak5PFL2L3')

### ntuple producer
process.ntupleProducer   = cms.EDAnalyzer('ntupleProducer',
  rootfilename      =    cms.untracked.string("nuTuple_TEST.root"),

  JetTag            =    cms.untracked.InputTag("ak5PFJets"),
  GenJetTag         =    cms.untracked.InputTag("ak5GenJets"),
  METTag            =    cms.untracked.InputTag("pfMet"),
  ElectronTag       =    cms.untracked.InputTag("gsfElectrons"),
  MuonTag           =    cms.untracked.InputTag("muons"),
  PhotonTag         =    cms.untracked.InputTag("photons"),
  TauTag            =    cms.untracked.InputTag("shrinkingConePFTauProducer"),
  PrimaryVtxTag     =    cms.untracked.InputTag("offlinePrimaryVerticesDAWithBS"),
  electronIDMap     =    cms.InputTag("simpleEleId80relIso"),
  rhoCorrTag	     =    cms.untracked.InputTag("kt6PFJets", "rho", "ntuples"),

  saveJets          =    cms.untracked.bool(True),
  saveElectrons     =    cms.untracked.bool(True),
  saveMuons         =    cms.untracked.bool(True),
  saveTaus          =    cms.untracked.bool(False),
  savePhotons       =    cms.untracked.bool(False),
  saveMET           =    cms.untracked.bool(True),
  saveGenJets       =    cms.untracked.bool(False),

  hltName           =    cms.untracked.string("HLT"),
  triggers          =    cms.untracked.vstring("HLT_Mu8_v",
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
                                               "HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v")
)

### Let it run
cmsSeq = cms.Sequence(
		process.offlinePrimaryVerticesDAWithBS 
      * process.simpleEleId80relIso
      * process.kt6PFJets
		* process.ak5PFJets
		* process.ak5PFJetsL1FastL2L3
		* process.ak5JetTracksAssociatorAtVertex 
		* process.btagging
      * process.HBHENoiseFilterResultProducer
		)

process.p = cms.Path(cmsSeq * process.ntupleProducer)

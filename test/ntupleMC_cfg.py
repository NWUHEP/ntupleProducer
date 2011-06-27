import FWCore.ParameterSet.Config as cms

process = cms.Process("ntuples") 

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'ERROR'

process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(100)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000))

process.load("Configuration.StandardSequences.Services_cff")
process.load('Configuration.StandardSequences.GeometryExtended_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

### Conditions tags
process.GlobalTag.globaltag = 'START42_V12::All' 

### Input files
#process.load("WJets")

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:/tmp/naodell/F88DBFB1-4B99-E011-AA88-0030486792BA.root'
#       '/store/data/Run2010A/JetMET/RECO/Sep17ReReco_v2/0026/FE9C387C-24C8-DF11-8CAA-003048679274.root'
)
)

# Select primary vertices
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

### ntuple producer
process.ntupleProducer   = cms.EDAnalyzer('ntupleProducer',

  RecoJetTag        =    cms.untracked.InputTag("ak5PFJets"),
  GenJetTag         =    cms.untracked.InputTag("ak5GenJets"),
  RecoMETTag        =    cms.untracked.InputTag("pfMet"),
  ElectronTag       =    cms.untracked.InputTag("gsfElectrons"),
  MuonTag           =    cms.untracked.InputTag("muons"),
  PrimaryVtxTag     =    cms.untracked.InputTag("offlinePrimaryVerticesDAWithBS"),
  gsfTrackTag       =    cms.untracked.InputTag("electronGsfTracks"),
  electronIDMap     =    cms.InputTag("simpleEleId80relIso"),
  trackTag          =    cms.untracked.InputTag("generalTracks"),
  barrelClusterTag  =    cms.untracked.InputTag("correctedHybridSuperClusters"),
  endcapClusterTag  =    cms.untracked.InputTag("correctedMulti5x5SuperClustersWithPreshower"),
  rhoCorrTag	     =    cms.untracked.InputTag("kt6PFJets", "rho", "ntuples"),
  nJets             =    cms.untracked.int32(2),
  nMuons            =    cms.untracked.int32(1),
  doPFJets          =    cms.untracked.bool(True),
  doGenJets         =    cms.untracked.bool(False),
  doFakeables       =    cms.untracked.bool(True),
  hltName           =    cms.untracked.string("HLT"),
  rootfilename      =    cms.untracked.string("nuTuple_WJets.root"),
  triggers          =    cms.untracked.vstring("HLT_Mu8_v", "HLT_Mu15_v", "HLT_Mu8_Jet40_v", "HLT_Mu13_Mu8_v", "HLT_Mu17_Mu8_v", "HLT_DoubleMu3_v", "HLT_DoubleMu6_v", "HLT_DoubleMu7_v", "HLT_Ele8_CaloIdL_CaloIsoVL_v", "HLT_Ele17_CaloIdL_CaloIsoVL_v", "HLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v", "HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v", "HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v")
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
		)

process.p = cms.Path(cmsSeq * process.ntupleProducer)


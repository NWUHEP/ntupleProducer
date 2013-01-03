import os
import FWCore.ParameterSet.Config as cms
from RecoEgamma.PhotonIdentification.isolationCalculator_cfi import *

process = cms.Process("PAT")

# real data or MC?
isRealData = False

# global tag
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')

if (isRealData):
    process.GlobalTag.globaltag = 'GR_R_53_V13::All'
else:
    process.GlobalTag.globaltag = 'START53_V11::All'

# Create good primary vertices for PF association
from PhysicsTools.SelectorUtils.pvSelector_cfi import pvSelector
process.goodOfflinePrimaryVertices = cms.EDFilter( "PrimaryVertexObjectFilter",
    filterParams = pvSelector.clone( minNdof = cms.double(4.0), maxZ = cms.double(24.0) ),
    src=cms.InputTag('offlinePrimaryVertices')
    )

# tau reconstruction configuration
#process.load("RecoTauTag.Configuration.RecoPFTauTag_cff")

# jet energy corrections
process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
process.load("CondCore.DBCommon.CondDBCommon_cfi")

# MET corrections
process.load('JetMETCorrections.Type1MET.pfMETCorrections_cff')
skipEMfractionThreshold = cms.double(0.90)
skipEM = cms.bool(True)
skipMuonSelection = cms.string('isGlobalMuon | isStandAloneMuon')
skipMuons = cms.bool(True)

if (isRealData):
    process.pfJetMETcorr.jetCorrLabel = cms.string("ak5PFL1FastL2L3Residual")

    # 53X b-jet discriminator calibration
    process.GlobalTag.toGet = cms.VPSet(
        cms.PSet(record = cms.string("BTagTrackProbability2DRcd"),
            tag = cms.string("TrackProbabilityCalibration_2D_Data53X_v2"),
            connect = cms.untracked.string("frontier://FrontierPrep/CMS_COND_BTAU")),
        cms.PSet(record = cms.string("BTagTrackProbability3DRcd"),
            tag = cms.string("TrackProbabilityCalibration_3D_Data53X_v2"),
            connect = cms.untracked.string("frontier://FrontierPrep/CMS_COND_BTAU"))
        )
else:
    # 53X b-jet discriminator calibration
    process.GlobalTag.toGet = cms.VPSet(
        cms.PSet(record = cms.string("BTagTrackProbability2DRcd"),
            tag = cms.string("TrackProbabilityCalibration_2D_Data53X_v2"),
            connect = cms.untracked.string("frontier://FrontierPrep/CMS_COND_BTAU")),
        cms.PSet(record = cms.string("BTagTrackProbability3DRcd"),
            tag = cms.string("TrackProbabilityCalibration_3D_Data53X_v2"),
            connect = cms.untracked.string("frontier://FrontierPrep/CMS_COND_BTAU"))
        )

### To get b-tags from ak5PFJets
process.load('RecoJets.JetAssociationProducers.ak5JTA_cff')
process.ak5JetTracksAssociatorAtVertex.jets = cms.InputTag("ak5PFJetsL1FastL2L3")
process.ak5JetTracksAssociatorAtCaloFace.jets = cms.InputTag("ak5PFJetsL1FastL2L3")
process.ak5JetExtender.jets = cms.InputTag("ak5PFJetsL1FastL2L3")

# jpt extras
process.load("RecoJets.Configuration.RecoPFJets_cff")
process.load("RecoJets.Configuration.RecoJPTJets_cff")
process.load("RecoJets.JetAssociationProducers.ak5JTA_cff")
process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
process.load('JetMETCorrections.Configuration.JetCorrectionServices_cff')


process.kt6PFJetsIso = process.kt6PFJets.clone()
process.kt6PFJetsIso.doRhoFastjet = True
process.kt6PFJetsIso.Rho_EtaMax = cms.double(2.5)


process.ak5JPTL1Offset.algorithm = 'AK5JPT'
process.ak5JetTracksAssociatorAtVertex.useAssigned = cms.bool(True)
process.ak5JetTracksAssociatorAtVertex.pvSrc = cms.InputTag("offlinePrimaryVertices")

process.jpt = cms.Sequence(
                        process.ak5JTA
                        * process.recoJPTJets
                        * process.ak5JPTJetsL1L2L3
                        * process.kt6PFJets
                        * process.ak5PFJetsL1FastL2L3
                        )

# pat sequences
process.load("PhysicsTools.PatAlgos.patSequences_cff")
from PhysicsTools.PatAlgos.patEventContent_cff import patEventContent
process.out = cms.OutputModule("PoolOutputModule",
                               fileName = cms.untracked.string('/tmp/patTuple.root'),
                               SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('ntuplePath')),
                               outputCommands = cms.untracked.vstring('keep *')#, *patEventContent )
                               )

from NWU.ntupleProducer.PatSequences_cff import addPatSequence
#addPatSequence(process, not isRealData, addPhotons = True)


# global options
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.cerr.FwkReport.reportEvery = 100

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
'''

process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(False),
                                     SkipEvent = cms.untracked.vstring('ProductNotFound')
                                    )

# event source
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
    #'/store/data/Run2012A/MuEG/AOD/13Jul2012-v1/0000/FEF59314-34D8-E111-8DF9-E0CB4E19F972.root'
    #'/store/mc/Summer12_DR53X/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S10_START53_V7A-v1/0002/D843FB2D-44D4-E111-A3C4-002481E75ED0.root'
    '/store/mc/Summer12_DR53X/GluGluToHToZG_M-125_8TeV-powheg-pythia6/AODSIM/PU_S10_START53_V7A-v1/0000/DEF04071-6EFA-E111-BA18-00266CFFC4D4.root'
    #'file:/tmp/naodell/TTJetsToHqToWWq_M-145_TuneZ2_8TeV_pythia6_v2_1_1_GPf.root'
)
)

# event counters
process.startCounter = cms.EDProducer("EventCountProducer")
process.endCounter = process.startCounter.clone()

#############################################
#### Met/Noise/BeamHalo filters  #############
## Following the recipe from this twiki:
## https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFilters
## Using *Tagging mode* for all filters (produces the boolean instead of filtering an event)
## Saving this boolean in the ntuples!
##############################################

process.load('CommonTools/RecoAlgos/HBHENoiseFilterResultProducer_cfi')
process.load("RecoMET.METFilters.hcalLaserEventFilter_cfi")
process.hcalLaserEventFilter.vetoByRunEventNumber=cms.untracked.bool(False)
process.hcalLaserEventFilter.vetoByHBHEOccupancy=cms.untracked.bool(True)
process.hcalLaserEventFilter.taggingMode = cms.bool(True)


process.load('RecoMET.METFilters.EcalDeadCellTriggerPrimitiveFilter_cfi')
## For AOD and RECO recommendation to use recovered rechits
process.EcalDeadCellTriggerPrimitiveFilter.tpDigiCollection = cms.InputTag("ecalTPSkimNA")
process.EcalDeadCellTriggerPrimitiveFilter.taggingMode = cms.bool(True)

# The section below is for the filter on Boundary Energy. Available in AOD in CMSSW>44x
# For releases earlier than 44x, one should make the following changes
# process.EcalDeadCellBoundaryEnergyFilter.recHitsEB = cms.InputTag("ecalRecHit","EcalRecHitsEB")
# process.EcalDeadCellBoundaryEnergyFilter.recHitsEE = cms.InputTag("ecalRecHit","EcalRecHitsEE")
process.load('RecoMET.METFilters.EcalDeadCellBoundaryEnergyFilter_cfi')
process.EcalDeadCellBoundaryEnergyFilter.taggingMode = cms.bool(True)
process.EcalDeadCellBoundaryEnergyFilter.cutBoundEnergyDeadCellsEB=cms.untracked.double(10)
process.EcalDeadCellBoundaryEnergyFilter.cutBoundEnergyDeadCellsEE=cms.untracked.double(10)
process.EcalDeadCellBoundaryEnergyFilter.cutBoundEnergyGapEB=cms.untracked.double(100)
process.EcalDeadCellBoundaryEnergyFilter.cutBoundEnergyGapEE=cms.untracked.double(100)
process.EcalDeadCellBoundaryEnergyFilter.enableGap=cms.untracked.bool(False)
process.EcalDeadCellBoundaryEnergyFilter.limitDeadCellToChannelStatusEB = cms.vint32(12,14)
process.EcalDeadCellBoundaryEnergyFilter.limitDeadCellToChannelStatusEE = cms.vint32(12,14)
# End of Boundary Energy filter configuration

# This one is not working for some reason:
#process.load('RecoMET.METFilters.trackingFailureFilter_cfi')
#process.trackingFailureFilter.taggingMode = cms.bool(True)

AllFilters = cms.Sequence(process.HBHENoiseFilterResultProducer
                          * process.hcalLaserEventFilter
                          * process.EcalDeadCellTriggerPrimitiveFilter
                          * process.EcalDeadCellBoundaryEnergyFilter
                          #* process.trackingFailureFilter
                            )


##### END OF Noise Filters ############

print '\n\nCommence ntuplization...\n\n'

### TFile service!
process.TFileService = cms.Service('TFileService',
                                  #fileName = cms.string('nuTuple.root')
                                   fileName = cms.string('nuTupleMVA.root')
                                   )

### ntuple producer
process.ntupleProducer   = cms.EDAnalyzer('ntupleProducer',

  printalot         =    cms.untracked.bool(False),

  photonIsoCalcTag  =    cms.PSet(isolationSumsCalculator),

  JetTag            =    cms.untracked.InputTag('ak5PFJetsL1FastL2L3'),
  GenJetTag         =    cms.untracked.InputTag('ak5GenJets'),
  #METTag            =    cms.untracked.InputTag('patMETsPFlow'),
  METTag            =    cms.untracked.InputTag('pfType1CorrectedMet'),
  ElectronTag       =    cms.untracked.InputTag('gsfElectrons'),
  MuonTag           =    cms.untracked.InputTag('muons'),
  PhotonTag         =    cms.untracked.InputTag('photons'),
  TauTag            =    cms.untracked.InputTag('selectedPatTausPFlow'),
  PrimaryVtxTag     =    cms.untracked.InputTag('offlinePrimaryVertices'),
  rhoCorrTag        =    cms.untracked.InputTag('kt6PFJets', 'rho', 'RECO'),
  rho25CorrTag      =    cms.untracked.InputTag('kt6PFJetsIso', 'rho', 'PAT'),
  rhoMuCorrTag      =    cms.untracked.InputTag('kt6PFJetsCentralNeutral', 'rho','RECO'),  # specifically for muon iso

  partFlowTag       =    cms.untracked.InputTag("particleFlow"), #,"Cleaned"),

  saveJets          =    cms.untracked.bool(False),
  saveElectrons     =    cms.untracked.bool(True),
  saveMuons         =    cms.untracked.bool(True),
  saveTaus          =    cms.untracked.bool(False),
  savePhotons       =    cms.untracked.bool(True),
  saveMET           =    cms.untracked.bool(True),
  saveGenJets       =    cms.untracked.bool(True),
  saveGenParticles  =    cms.untracked.bool(True),


  ecalTPFilterTag    =    cms.untracked.InputTag("EcalDeadCellTriggerPrimitiveFilter",""),
  ecalBEFilterTag    =    cms.untracked.InputTag("EcalDeadCellBoundaryEnergyFilter",""),
  hcalHBHEFilterTag  =    cms.untracked.InputTag("HBHENoiseFilterResultProducer","HBHENoiseFilterResult"),
  hcalLaserFilterTag =    cms.untracked.InputTag("hcalLaserEventFilter",""),

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
                                               "HLT_Mu17_TkMu8_v",
                                               "HLT_Mu22_TkMu8_v",
                                               "HLT_Mu22_TkMu22_v",

                                               "HLT_Ele8_CaloIdL_CaloIsoVL_v",
                                               "HLT_Ele17_CaloIdL_CaloIsoVL_v",
                                               "HLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v",
                                               "HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v",
                                               "HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v",
                                               "HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v",
                                               "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v",

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
                                               "HLT_Photon150_v",
                                               "HLT_Photon22_R9Id90_HE10_Iso40_EBOnly_v",
                                               "HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_v",
                                               "HLT_Photon50_R9Id90_HE10_Iso40_EBOnly_v",
                                               "HLT_Photon75_R9Id90_HE10_Iso40_EBOnly_v",
                                               "HLT_Photon90_R9Id90_HE10_Iso40_EBOnly_v",

                                               "HLT_Mu17_Ele8_CaloIdL_v",
                                               "HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v",
                                               "HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v",
                                               "HLT_Mu3_Ele8_CaloIdT_CaloIsoVL_v",
                                               "HLT_Mu8_Ele17_CaloIdL_v",
                                               "HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v",
                                               "HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v",

                                               "HLT_Ele18_CaloIdVT_TrkIdT_MediumIsoPFTau20_v",
                                               "HLT_Ele20_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_MediumIsoPFTau20_v",
                                               "HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_MediumIsoPFTau25_v",
                                               "HLT_HT400_DoubleIsoPFTau10_Trk3_PFMHT50_v",
                                               "HLT_IsoMu15_eta2p1_LooseIsoPFTau20_v",
                                               "HLT_IsoMu15_eta2p1_MediumIsoPFTau20_v",
                                               "HLT_IsoMu15_eta2p1_TightIsoPFTau20_v",
                                               "HLT_Mu15_LooseIsoPFTau15_v"
)
)

process.ntuplePath = cms.Path(
        process.goodOfflinePrimaryVertices
        * process.producePFMETCorrections
        #* process.PFTau
        #* process.patDefaultSequence
        #* process.jpt
        * process.kt6PFJetsIso
        * process.ak5PFJetsL1FastL2L3
        * process.ak5JetTracksAssociatorAtVertex
        * process.btagging
        * AllFilters
        * process.ntupleProducer
        )

#process.outpath = cms.EndPath(process.out)

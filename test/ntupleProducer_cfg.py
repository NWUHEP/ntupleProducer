import os
import FWCore.ParameterSet.Config as cms
from RecoEgamma.PhotonIdentification.isolationCalculator_cfi import *

process = cms.Process("NTUPLE")

# real data or MC?
isRealData = False

# global tag
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')

if (isRealData):
    process.GlobalTag.globaltag = 'GR_P_V42_AN4::All'
else:
    process.GlobalTag.globaltag = 'START53_V19D::All'

# Create good primary vertices for PF association
from PhysicsTools.SelectorUtils.pvSelector_cfi import pvSelector
process.goodOfflinePrimaryVertices = cms.EDFilter( "PrimaryVertexObjectFilter",
    filterParams = pvSelector.clone( minNdof = cms.double(4.0), maxZ = cms.double(24.0) ),
    src=cms.InputTag('offlinePrimaryVertices')
    )

# jet energy corrections
process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
process.load("CondCore.DBCommon.CondDBCommon_cfi")


# Add MET collection for PAT
#from PhysicsTools.PatAlgos.tools.metTools import *
#addPfMET(process,'PF')
#addTcMET(process,"TC")

# MET corrections Type 1 and x,y corrections
process.load('JetMETCorrections.Type1MET.pfMETCorrections_cff')
process.load("JetMETCorrections.Type1MET.pfMETsysShiftCorrections_cfi")

# use for 2012 Data
if (isRealData):
    process.pfMEtSysShiftCorr.parameter = process.pfMEtSysShiftCorrParameters_2012runAvsNvtx_data
# use for Spring'12 MC
else:
    process.pfMEtSysShiftCorr.parameter = process.pfMEtSysShiftCorrParameters_2012runAvsNvtx_mc

process.pfType1CorrectedMet.srcType1Corrections = cms.VInputTag(
    cms.InputTag('pfJetMETcorr', 'type1') ,
    cms.InputTag('pfMEtSysShiftCorr')
)
process.pfType1p2CorrectedMet.srcType1Corrections = cms.VInputTag(
    cms.InputTag('pfJetMETcorr', 'type1') ,
    cms.InputTag('pfMEtSysShiftCorr')
)

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

## The iso-based HBHE noise filter
process.load('CommonTools/RecoAlgos/HBHENoiseFilterResultProducer_cfi')

## This might need updating for the VBF Parked - waiting for the recommendation
process.load("RecoMET.METFilters.hcalLaserEventFilter_cfi")
process.hcalLaserEventFilter.vetoByRunEventNumber=cms.untracked.bool(False)
process.hcalLaserEventFilter.vetoByHBHEOccupancy=cms.untracked.bool(True)
process.hcalLaserEventFilter.taggingMode = cms.bool(True)

## Ecal Dead Cell Filter
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


## The Good vertices collection needed by the tracking failure filter
process.goodVertices = cms.EDFilter(
      "VertexSelector",
        filter = cms.bool(False),
        src = cms.InputTag("offlinePrimaryVertices"),
        cut = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && position.rho < 2")
      )
## The tracking failure filter
process.load('RecoMET.METFilters.trackingFailureFilter_cfi')
process.trackingFailureFilter.taggingMode = cms.bool(True)


#Bad EE SC filter, not needed but goot to have them
process.load('RecoMET.METFilters.eeBadScFilter_cfi')
process.eeBadScFilter.tagginMode = cms.bool(True)

## The tracking POG filters
process.load('RecoMET.METFilters.trackingPOGFilters_cff')
## NOTE: to make tagging mode of the tracking POG filters (three of them), please do:
process.manystripclus53X.taggedMode = cms.untracked.bool(True)
process.manystripclus53X.forcedValue = cms.untracked.bool(False)
process.toomanystripclus53X.taggedMode = cms.untracked.bool(True)
process.toomanystripclus53X.forcedValue = cms.untracked.bool(False)
process.logErrorTooManyClusters.taggedMode = cms.untracked.bool(True)
process.logErrorTooManyClusters.forcedValue = cms.untracked.bool(False)
## Also the stored boolean for the three filters is opposite to what we usually
## have for other filters, i.e., true means rejected bad events while false means
## good events.


AllFilters = cms.Sequence(process.HBHENoiseFilterResultProducer
                          * process.hcalLaserEventFilter
                          * process.EcalDeadCellTriggerPrimitiveFilter
                          * process.EcalDeadCellBoundaryEnergyFilter
                          * process.goodVertices * process.trackingFailureFilter
                          * process.eeBadScFilter
                          #* process.trkPOGFilters
                          * ~process.manystripclus53X #trkPOGFilter1
                          * ~process.toomanystripclus53X #trkPOGFilter2
                          * ~process.logErrorTooManyClusters #trkPOGFilter 3
                          )

# event source
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(10000))
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
      '/store/user/bpollack/POWHEG_PYTHIA8_H_Zg_8TeV_cff_py_GEN_SIM_test/POWHEG_PYTHIA8_H_Zg_8TeV_cff_py_GEN_SIM_REDIGI_DIGI_L1_DIGI2RAW_HLT_PU_RAW2DIGI_L1Reco_RECO_VALIDATION_DQM_PU_test/fc2433b7e79d5e0440261850868ad62f/POWHEG_PYTHIA8_H_Zg_8TeV_cff_py_GEN_SIM_REDIGI_DIGI_L1_DIGI2RAW_HLT_PU_RAW2DIGI_L1Reco_RECO_VALIDATION_DQM_PU_test_1_1_cJt.root'
      #'file:/uscms_data/d2/bpollack/genProd/CMSSW_5_3_8/src/test/STEP2_RAW2DIGI_L1Reco_RECO_VALIDATION_DQM_PU.root'
      #'/store/data/Run2012A/SingleMu/AOD/13Jul2012-v1/00000/82CD2186-C7D0-E111-8F3C-0017A4771018.root'
      #'/store/mc/Summer12_DR53X/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S10_START53_V7A-v1/0002/D843FB2D-44D4-E111-A3C4-002481E75ED0.root'
    #'/store/data/Run2012D/DoubleMu/AOD/PromptReco-v1/000/208/341/285B355D-553D-E211-A3FC-BCAEC532971E.root'
    #'file:/tmp/naodell/TTJetsToHqToWWq_M-125_TuneZ2_8TeV_pythia6_v2_1_1_p64.root'
)
)


##### END OF Noise Filters ############

print '\n\nCommence ntuplization...\n\n'

### TFile service!
process.TFileService = cms.Service('TFileService',
                                  fileName = cms.string('nuTuple.root')
                                   )

### pfNoPU Sequence for electron MVA
process.pfPileUp = cms.EDProducer("PFPileUp",
    PFCandidates = cms.InputTag("particleFlow"),
    Enable = cms.bool(True),
    checkClosestZVertex = cms.bool(True),
    verbose = cms.untracked.bool(False),
    Vertices = cms.InputTag("offlinePrimaryVertices")
)

process.pfNoPileUp = cms.EDProducer("TPPFCandidatesOnPFCandidates",
    bottomCollection = cms.InputTag("particleFlow"),
    enable = cms.bool(True),
    topCollection = cms.InputTag("pfPileUp"),
    name = cms.untracked.string('pileUpOnPFCandidates'),
    verbose = cms.untracked.bool(False)
)

process.pfNoPUSeq = cms.Sequence(process.pfPileUp + process.pfNoPileUp)


### ntuple producer
process.ntupleProducer   = cms.EDAnalyzer('ntupleProducer',

  verboseTrigs         =    cms.untracked.bool(False),
  verboseMVAs          =    cms.untracked.bool(False),

  photonIsoCalcTag  =    cms.PSet(isolationSumsCalculator),

  JetTag            =    cms.untracked.InputTag('ak5PFJetsL1FastL2L3'),
  GenJetTag         =    cms.untracked.InputTag('ak5GenJets'),
  METTag            =    cms.untracked.InputTag('pfType1CorrectedMet'),
  ElectronTag       =    cms.untracked.InputTag('gsfElectrons'),
  MuonTag           =    cms.untracked.InputTag('muons'),
  PhotonTag         =    cms.untracked.InputTag('photons'),
  PrimaryVtxTag     =    cms.untracked.InputTag('offlinePrimaryVertices'),
  rhoCorrTag        =    cms.untracked.InputTag('kt6PFJets', 'rho', 'RECO'),
  rho25CorrTag      =    cms.untracked.InputTag('kt6PFJetsIso', 'rho', 'NTUPLE'),
  rhoMuCorrTag      =    cms.untracked.InputTag('kt6PFJetsCentralNeutral', 'rho','RECO'),  # specifically for muon iso

  partFlowTag       =    cms.untracked.InputTag("particleFlow"), #,"Cleaned"),

  saveJets          =    cms.untracked.bool(True),
  saveElectrons     =    cms.untracked.bool(True),
  saveMuons         =    cms.untracked.bool(True),
  savePhotons       =    cms.untracked.bool(True),
  saveMET           =    cms.untracked.bool(True),
  saveGenJets       =    cms.untracked.bool(True),
  saveGenParticles  =    cms.untracked.bool(True),

  ecalTPFilterTag    =    cms.untracked.InputTag("EcalDeadCellTriggerPrimitiveFilter",""),
  ecalBEFilterTag    =    cms.untracked.InputTag("EcalDeadCellBoundaryEnergyFilter",""),
  hcalHBHEFilterTag  =    cms.untracked.InputTag("HBHENoiseFilterResultProducer","HBHENoiseFilterResult"),
  hcalLaserFilterTag =    cms.untracked.InputTag("hcalLaserEventFilter",""),
  trackingFailureTag =    cms.untracked.InputTag("trackingFailureFilter",""),
  eeBadScFilterTag   =    cms.untracked.InputTag("eeBadScFilter",""),
  trkPOGFiltersTag1  =    cms.untracked.InputTag("manystripclus53X",""),
  trkPOGFiltersTag2  =    cms.untracked.InputTag("toomanystripclus53X",""),
  trkPOGFiltersTag3  =    cms.untracked.InputTag("logErrorTooManyClusters",""),

  hltName           =    cms.untracked.string("HLT"),
  triggers          =    cms.untracked.vstring(
                                               "HLT_IsoMu24_v",
                                               "HLT_IsoMu24_eta2p1_v",
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

                                               "HLT_DiJet20_MJJ650_AllJets_DEta3p5_HT120_VBF_v"
                                               "HLT_DiJet30_MJJ700_AllJets_DEta3p5_VBF_v"
                                               "HLT_DiJet35_MJJ650_AllJets_DEta3p5_VBF_v"
                                               "HLT_DiJet35_MJJ700_AllJets_DEta3p5_VBF_v"
                                               "HLT_DiJet35_MJJ750_AllJets_DEta3p5_VBF_v"

                                               "HLT_Photon30_R9Id90_CaloId_HE10_Iso40_EBOnly_Met25_HBHENoiseCleaned_v"
                                               "HLT_Photon30_R9Id90_CaloId_HE10_Iso40_EBOnly_v"
                                               "HLT_Photon30_v"
)
)

process.ntuplePath = cms.Path(
        process.goodOfflinePrimaryVertices
        * process.pfMEtSysShiftCorrSequence
        * process.producePFMETCorrections
        * process.pfNoPUSeq
        #* process.patDefaultSequence
        * process.kt6PFJetsIso
        * process.ak5PFJetsL1FastL2L3
        * process.ak5JetTracksAssociatorAtVertex
        * process.btagging
        * AllFilters
        * process.ntupleProducer
        )

#process.outpath = cms.EndPath(process.out)

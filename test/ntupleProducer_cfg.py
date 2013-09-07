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
    process.GlobalTag.globaltag = 'FT_53_V21_AN3::All'
    process.load('JetMETCorrections.METPUSubtraction.mvaPFMET_leptons_data_cff')
else:
    process.GlobalTag.globaltag = 'START53_V27::All'
    process.load('JetMETCorrections.METPUSubtraction.mvaPFMET_leptons_cff')

# Create good primary vertices for PF association
from PhysicsTools.SelectorUtils.pvSelector_cfi import pvSelector
process.goodOfflinePrimaryVertices = cms.EDFilter( "PrimaryVertexObjectFilter",
    filterParams = pvSelector.clone( minNdof = cms.double(4.0), maxZ = cms.double(24.0) ),
    src=cms.InputTag('offlinePrimaryVertices') #Standard Primary Vertex Collection
#    src=cms.InputTag('offlinePrimaryVerticesWithBS') #Primary Vertices Collection constrained by beamspot
    )

# jet energy corrections
process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
process.load("CondCore.DBCommon.CondDBCommon_cfi")


# Add MET collection for PAT
#from PhysicsTools.PatAlgos.tools.metTools import *
#addPfMET(process,'PF')
#addTcMET(process,"TC")

##Testing MET significance
process.load("RecoMET.METProducers.PFMET_cfi")
process.pfMet1 = process.pfMet.clone(alias="PFMET1")


# MET corrections Type 1 and x,y corrections
process.load('JetMETCorrections.Type1MET.pfMETCorrections_cff')
process.load("JetMETCorrections.Type1MET.pfMETsysShiftCorrections_cfi")
process.load("JetMETCorrections.Type1MET.pfMETCorrectionType0_cfi")


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

process.pfType1CorrectedMetType0.srcType1Corrections = cms.VInputTag(
    cms.InputTag('pfMETcorrType0'),
    cms.InputTag('pfJetMETcorr', 'type1') ,
    cms.InputTag('pfMEtSysShiftCorr')
)

process.load("RecoMET.METProducers.pfChargedMET_cfi")
process.load("RecoMET.METProducers.TrackMET_cfi")

#if (isRealData == False):
#process.load("RecoMET.Configuration.GenMETParticles_cff")
#process.load("RecoMET.METProducers.genMetCalo_cfi")
#process.load("RecoMET.METProducers.MetMuonCorrections_cff")

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
#process.load("RecoJets.Configuration.RecoPFJets_cff")
#process.load("RecoJets.Configuration.RecoJPTJets_cff")
#process.load("RecoJets.JetAssociationProducers.ak5JTA_cff")
#process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
#process.load('JetMETCorrections.Configuration.JetCorrectionServices_cff')


process.kt6PFJetsIso = process.kt6PFJets.clone()
process.kt6PFJetsIso.doRhoFastjet = True
process.kt6PFJetsIso.Rho_EtaMax = cms.double(2.5)


#process.ak5JPTL1Offset.algorithm = 'AK5JPT'
#process.ak5JetTracksAssociatorAtVertex.useAssigned = cms.bool(True)
#process.ak5JetTracksAssociatorAtVertex.pvSrc = cms.InputTag("offlinePrimaryVertices")
#process.ak5JetTracksAssociatorAtVertex.pvSrc = cms.InputTag("offlinePrimaryVerticesWithBS")

#process.jpt = cms.Sequence(
#                        process.ak5JTA
#                        * process.recoJPTJets
#                        * process.ak5JPTJetsL1L2L3
#                        * process.kt6PFJets
#                        * process.ak5PFJetsL1FastL2L3
#                        )

# for jet pileup ID variables
from RecoJets.JetProducers.PileupJetIDParams_cfi import *

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
#        src = cms.InputTag("offlinePrimaryVerticesWithBS"),
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
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(100))
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
	#'/store/mc/Summer12_DR53X/GluGluToHToWWTo2LAndTau2Nu_M-125_8TeV-powheg-pythia6/AODSIM/PU_S10_START53_V7A-v1/0000/DE5F727F-8BFC-E111-8576-002618FDA263.root'
#	'root://eoscms//eos/cms/store/user/cmkuo/GluGluToHToZG_M-125_8TeV-powheg-pythia6/HZg_nunug_ggH_m125_RECO_v1/3664d28163503ca8171ba37083c39fc9/STEP2_RAW2DIGI_L1Reco_RECO_PU_100_1_fXq.root'
#'/store/data/Run2012D/SinglePhotonParked/AOD/22Jan2013-v1/30004/144D7268-4086-E211-9DC1-001E673984C1.root'
#'/store/mc/Summer12_DR53X/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S10_START53_V7A-v1/0002/D843FB2D-44D4-E111-A3C4-002481E75ED0.root'
    #'/store/data/Run2012D/DoubleMu/AOD/PromptReco-v1/000/208/341/285B355D-553D-E211-A3FC-BCAEC532971E.root'
    #'file:/tmp/naodell/TTJetsToHqToWWq_M-125_TuneZ2_8TeV_pythia6_v2_1_1_p64.root'
    'file:/uscms_data/d2/bpollack/genProd/CMSSW_5_3_8/src/test/testOut2_v2/PYTHIA8_175_POWHEG_H_Zg_8TeV_cff_py_GEN_SIM_REDIGI_DIGI_L1_DIGI2RAW_HLT_PU_STEP2_RAW2DIGI_L1Reco_RECO_VALIDATION_DQM_PU_50.root'
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
#    Vertices = cms.InputTag("offlinePrimaryVerticesWithBS")
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
  jetPUIdAlgo       =    cms.PSet(full_5x),

  JetTag            =    cms.untracked.InputTag('ak5PFJetsL1FastL2L3'),
  JecTag            =    cms.string("AK5PF"),
  GenJetTag         =    cms.untracked.InputTag('ak5GenJets'),
  METTag            =    cms.untracked.InputTag('pfType1CorrectedMet'),
  TrackMETTag       =    cms.untracked.InputTag('trackMet'),
  ElectronTag       =    cms.untracked.InputTag('gsfElectrons'),
  MuonTag           =    cms.untracked.InputTag('muons'),
  PhotonTag         =    cms.untracked.InputTag('photons'),
  PrimaryVtxTag     =    cms.untracked.InputTag('offlinePrimaryVertices'),
  rhoCorrTag        =    cms.untracked.InputTag('kt6PFJets', 'rho', 'RECO'),
  rho25CorrTag      =    cms.untracked.InputTag('kt6PFJetsIso', 'rho', 'NTUPLE'),
  rhoMuCorrTag      =    cms.untracked.InputTag('kt6PFJetsCentralNeutral', 'rho','RECO'),  # specifically for muon iso

  T0METTag	    =	 cms.untracked.InputTag('pfType1CorrectedMetType0'),
  T2METTag	    =	 cms.untracked.InputTag('pfType1p2CorrectedMet'),

  partFlowTag       =    cms.untracked.InputTag("particleFlow"), #,"Cleaned"),

  saveJets          =    cms.untracked.bool(True),
  saveElectrons     =    cms.untracked.bool(True),
  saveMuons         =    cms.untracked.bool(True),
  savePhotons       =    cms.untracked.bool(True),
  saveMET           =    cms.untracked.bool(True),
  saveGenJets       =    cms.untracked.bool(True),
  saveGenParticles  =    cms.untracked.bool(True),
  saveTrackMET      =    cms.untracked.bool(True),
  saveT0MET	    =    cms.untracked.bool(True),
  saveT2MET	    =    cms.untracked.bool(True),

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
                                               "HLT_Mu13_Mu8_v",
                                               "HLT_Mu17_Mu8_v",
                                               "HLT_DoubleMu7_v",
                                               "HLT_Mu17_TkMu8_v",
                                               "HLT_Mu22_TkMu8_v",
                                               "HLT_Mu22_TkMu22_v",

                                               "HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v",
                                               "HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v",
                                               "HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v",
                                               "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v",

                                               "HLT_Mu17_Ele8_CaloIdL_v",
                                               "HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v",
                                               "HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v",
                                               "HLT_Mu8_Ele17_CaloIdL_v",
                                               "HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v",
                                               "HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v",


                                               "HLT_Photon30_R9Id90_CaloId_HE10_Iso40_EBOnly_Met25_HBHENoiseCleaned",
                                               "HLT_Photon30_R9Id90_CaloId_HE10_Iso40_EBOnly",
                                               "HLT_Photon30",
                                               "HLT_DiJet20_MJJ650_AllJets_DEta3p5_HT120_VBF",
                                               "HLT_DiJet30_MJJ700_AllJets_DEta3p5_VBF",
                                               "HLT_DiJet35_MJJ650_AllJets_DEta3p5_VBF",
                                               "HLT_DiJet35_MJJ700_AllJets_DEta3p5_VBF",
                                               "HLT_DiJet35_MJJ750_AllJets_DEta3p5_VBF"
)
)

process.ntuplePath = cms.Path(
      process.goodOfflinePrimaryVertices
		* process.type0PFMEtCorrection
    * process.pfMEtSysShiftCorrSequence
    * process.producePFMETCorrections
    * process.pfNoPUSeq
		* process.particleFlowForChargedMET
		* process.pfChargedMET
		* process.trackMet
   #* process.patDefaultSequence
    * process.kt6PFJetsIso
    * process.ak5PFJetsL1FastL2L3
    * process.ak5JetTracksAssociatorAtVertex
    * process.btagging
    * AllFilters
    * process.pfMEtMVAsequence
    * process.pfMet1
    * process.ntupleProducer
)

#if (isRealData == False):
#	process.ntuplePath = cms.Path(
#        	process.goodOfflinePrimaryVertices
#		* process.genCandidatesForMET
#		* process.corMetGlobalMuons
#		* process.genMetCalo
#	        * process.pfMEtSysShiftCorrSequence
#	        * process.producePFMETCorrections
#	        * process.pfNoPUSeq
#		* process.particleFlowForChargedMET
#		* process.pfChargedMET
#		* process.trackMet
#	        #* process.patDefaultSequence
#	        * process.kt6PFJetsIso
#	        * process.ak5PFJetsL1FastL2L3
#	        * process.ak5JetTracksAssociatorAtVertex
#	        * process.btagging
#	        * AllFilters
#	        * process.ntupleProducer
#	        )
#process.outpath = cms.EndPath(process.out)

import os
import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
from RecoEgamma.PhotonIdentification.isolationCalculator_cfi import *

process = cms.Process("NTUPLE")

options = VarParsing.VarParsing ('analysis')
options.maxEvents = 10
#options.inputFiles= '/store/data/Run2012C/SingleMu/AOD/22Jan2013-v1/30010/C0E05558-9078-E211-9E02-485B39800B65.root'
#options.inputFiles= '/store/data/Run2012C/DoublePhoton/AOD/22Jan2013-v2/30001/72DE4526-F370-E211-B370-00304867920A.root'
#options.loadFromFile('inputFiles','PYTHIA8_175_H_Zg_8TeV.txt')
#options.loadFromFile('inputFiles','PYTHIA8_175_POWHEG_PDF7_H_Zg_8TeV.txt')
#options.inputFiles = '/store/data/Run2012A/DoubleElectron/AOD/13Jul2012-v1/00000/00347915-EED9-E111-945A-848F69FD2817.root'
#options.inputFiles = '/store/user/andrey/Higgs_To_MuMuGamma_Dalitz_MH125_Mll_0to50_MadgraphHEFT_pythia6/AODSIM_v2/39bf61f738ba3bdb8860f0848073cc88/aodsim_100_1_hLi.root'
#'/store/mc/Summer12_DR53X/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S10_START53_V7A-v1/0002/D843FB2D-44D4-E111-A3C4-002481E75ED0.root'
#options.inputFiles = 'file:/uscms_data/d2/bpollack/genProd/CMSSW_5_3_8/src/test/testOut2_v2/PYTHIA8_175_POWHEG_H_Zg_8TeV_cff_py_GEN_SIM_REDIGI_DIGI_L1_DIGI2RAW_HLT_PU_STEP2_RAW2DIGI_L1Reco_RECO_VALIDATION_DQM_PU_50.root'
options.inputFiles = '/store/mc/Summer12_DR53X/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S10_START53_V7A-v1/0000/02CDCF05-BED2-E111-85F4-0030486740BA.root'
#options.inputFiles = '/store/user/andrey/MCFM_lord_hzgamma_8TeV_LHE_pythia6_v2/AODSIM/39bf61f738ba3bdb8860f0848073cc88/aodsim_100_1_BGG.root'
#options.inputFiles = '/store/data/Run2012D/SinglePhotonParked/AOD/22Jan2013-v1/30004/144D7268-4086-E211-9DC1-001E673984C1.root'
#options.inputFiles = 'file:/uscms_data/d2/bpollack/genProd/CMSSW_5_3_8/src/test/testOut2_v2/PYTHIA8_175_POWHEG_H_Zg_8TeV_cff_py_GEN_SIM_REDIGI_DIGI_L1_DIGI2RAW_HLT_PU_STEP2_RAW2DIGI_L1Reco_RECO_VALIDATION_DQM_PU_8.root'

options.register("isRealData",
                 0,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 "0 if running on MC and 1 if running on Data")

options.parseArguments()

# real data or MC?
isRealData = options.isRealData

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

#commenting out  pat sequences. we don't use them
#process.load("PhysicsTools.PatAlgos.patSequences_cff")
#from PhysicsTools.PatAlgos.patEventContent_cff import patEventContent

#from NWU.ntupleProducer.PatSequences_cff import addPatSequence
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
##### END OF Noise Filters ############


# Electron MVA ID producer:
process.load('EgammaAnalysis/ElectronTools/electronIdMVAProducer_cfi')

# Electron Regression (post moriond recommendation)
process.load('EgammaAnalysis/ElectronTools/electronRegressionEnergyProducer_cfi')
process.eleRegressionEnergy.inputElectronsTag    = cms.InputTag('gsfElectrons')
process.eleRegressionEnergy.inputCollectionType  = cms.uint32(0)
process.eleRegressionEnergy.useRecHitCollections = cms.bool(True)
process.eleRegressionEnergy.produceValueMaps     = cms.bool(True)
process.eleRegressionEnergy.energyRegressionType = cms.uint32(2)
process.eleRegressionEnergy.regressionInputFile = cms.string("EgammaAnalysis/ElectronTools/data/eleEnergyRegWeights_WithSubClusters_VApr15.root")

# Electron Combination for calibration (post moriond recommendation)
process.load('EgammaAnalysis/ElectronTools/calibratedElectrons_cfi')
process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
    calibratedElectrons = cms.PSet(
      initialSeed = cms.untracked.uint32(1),
      engineName = cms.untracked.string('TRandom3')
    ),
)
if (isRealData):
  process.calibratedElectrons.isMC = cms.bool(False)
  process.calibratedElectrons.inputDataset = cms.string("22Jan2013ReReco")
else:
  process.calibratedElectrons.isMC = cms.bool(True)
  process.calibratedElectrons.inputDataset = cms.string("Summer12_LegacyPaper")

process.calibratedElectrons.updateEnergyError = cms.bool(True)
process.calibratedElectrons.correctionsType   = cms.int32(2)
process.calibratedElectrons.combinationType   = cms.int32(3)
process.calibratedElectrons.lumiRatio         = cms.double(1.0)
process.calibratedElectrons.verbose           = cms.bool(False)
process.calibratedElectrons.synchronization   = cms.bool(False)
process.calibratedElectrons.applyLinearityCorrection = cms.bool(True)
#process.calibratedElectrons.scaleCorrectionsInputPath = cms.string("EgammaAnalysis/ElectronTools/data/scalesMoriond.csv")
#process.calibratedElectrons.combinationRegressionInputPath = cms.string("EgammaAnalysis/ElectronTools/data/eleEnergyReg2012Weights_V1.root")


# event source
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(options.maxEvents))
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(options.inputFiles)
)


print '\n\nCommence ntuplization...\n\n'

### TFile service!
process.TFileService = cms.Service('TFileService',
    #fileName = cms.string('nuTuple.root')
    fileName = cms.string('~/EOS/V08_01_8TeV/ggHZG_M125_Pythia8_175_POWHEG_PDF7/nuTuple_9.root')
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



from SHarper.HEEPAnalyzer.HEEPSelectionCuts_cfi import *
process.heepIdNoIso = cms.EDProducer("HEEPIdValueMapProducer",
                                     eleLabel = cms.InputTag("gsfElectrons"),
                                     barrelCuts = cms.PSet(heepBarrelCuts),
                                     endcapCuts = cms.PSet(heepEndcapCuts),
                                     eleIsolEffectiveAreas = cms.PSet(heepEffectiveAreas),
                                     eleRhoCorrLabel = cms.InputTag("kt6PFJets", "rho"),
                                     verticesLabel = cms.InputTag("offlinePrimaryVertices"),
                                     applyRhoCorrToEleIsol = cms.bool(True),
                                     writeIdAsInt = cms.bool(True)
                                     )
process.heepIdNoIso.barrelCuts.cuts=cms.string("et:detEta:ecalDriven:dEtaIn:dPhiIn:hadem:e2x5Over5x5:nrMissHits:dxy")
process.heepIdNoIso.endcapCuts.cuts=cms.string("et:detEta:ecalDriven:dEtaIn:dPhiIn:hadem:sigmaIEtaIEta:nrMissHits:dxy")

process.heepIdNoIsoEles = cms.EDProducer("tsw::HEEPGsfProducer", cutValueMap = cms.InputTag("heepIdNoIso"),
                                         inputGsfEles = cms.InputTag("gsfElectrons")  )

# Boosted Z ModEleIso: 1b) Calculating the modified iso. values using BstdZeeTools EDProducer

from TSWilliams.BstdZeeTools.bstdzeemodisolproducer_cff import *
process.modElectronIso = cms.EDProducer("BstdZeeModIsolProducer",
                                              bstdZeeModIsolParams, vetoGsfEles = cms.InputTag("heepIdNoIsoEles") )


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
  rhoCorrTag        =    cms.untracked.InputTag('kt6PFJets',    'rho', 'RECO'),
  rho25CorrTag      =    cms.untracked.InputTag('kt6PFJetsIso', 'rho', 'NTUPLE'),
  rhoMuCorrTag      =    cms.untracked.InputTag('kt6PFJetsCentralNeutral', 'rho','RECO'),  # specifically for muon iso

  T0METTag	    =	 cms.untracked.InputTag('pfType1CorrectedMetType0'),
  T2METTag	    =	 cms.untracked.InputTag('pfType1p2CorrectedMet'),

  partFlowTag       =  cms.untracked.InputTag("particleFlow"), #,"Cleaned"),
  skimLepton        =  cms.untracked.bool(False),

  saveMuons         =    cms.untracked.bool(True),
  saveJets          =    cms.untracked.bool(True),
  saveElectrons     =    cms.untracked.bool(True),
  saveEleCrystals   =    cms.untracked.bool(True),
  savePhotons       =    cms.untracked.bool(True),
  savePhoCrystals   =    cms.untracked.bool(True),
  saveMoreEgammaVars=    cms.untracked.bool(True),

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

  #for SC footprint removal

  isolation_cone_size_forSCremoval = cms.untracked.double(0.3),


  hltName           =    cms.untracked.string("HLT"),
  triggers          =    cms.untracked.vstring(
                                               "HLT_Mu13_Mu8_v",
                                               "HLT_Mu17_Mu8_v",
                                               "HLT_Mu17_TkMu8_v",
                                               "HLT_Mu22_TkMu8_v",
                                               "HLT_Mu22_TkMu22_v",
                                               "HLT_IsoMu24_v",
                                               "HLT_IsoMu24_eta2p1_v",

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

                                               "HLT_Mu22_Photon22_CaloIdL_v",
                                               "HLT_Ele27_WP80_v",
                                               "HLT_Photon36_CaloId10_Iso50_Photon22_CaloId10_Iso50_v",
                                               "HLT_Photon36_CaloId10_Iso50_Photon22_R9Id85_v",
                                               "HLT_Photon36_R9Id85_OR_CaloId10_Iso50_Photon22_R9Id85_OR_CaloId10_Iso50_v",
                                               "HLT_Photon36_R9Id85_Photon22_CaloId10_Iso50_v",
                                               "HLT_Photon36_R9Id85_Photon22_R9Id85_v",


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

    * process.eleRegressionEnergy
    * process.calibratedElectrons
    * process.mvaTrigV0

    * process.heepIdNoIso
    * process.heepIdNoIsoEles
    * process.modElectronIso

    * process.ntupleProducer
)

#process.out = cms.OutputModule("PoolOutputModule",
#                               fileName = cms.untracked.string('/tmp/myTuple.root'),
#                               SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('ntuplePath')),
#                               outputCommands = cms.untracked.vstring('keep *')
#                               )
#process.outpath = cms.EndPath(process.out)

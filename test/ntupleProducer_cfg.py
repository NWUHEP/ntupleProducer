import os
import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
from RecoEgamma.PhotonIdentification.isolationCalculator_cfi import *

process = cms.Process("NTUPLE")

options = VarParsing.VarParsing ('analysis')
options.maxEvents = 100
options.inputFiles = '/store/mc/Summer12_DR53X/TTH_HToZG_M-135_8TeV-pythia8175/AODSIM/PU_RD1_START53_V7N-v2/00000/0A4C1013-9287-E311-937C-003048D4397E.root'
#options.inputFiles= '/store/data/Run2012C/SingleMu/AOD/22Jan2013-v1/30010/C0E05558-9078-E211-9E02-485B39800B65.root'
#options.loadFromFile('inputFiles','PYTHIA8_175_H_Zg_8TeV.txt')
#options.loadFromFile('inputFiles','PYTHIA8_175_POWHEG_PDF7_H_Zg_8TeV.txt')
#options.inputFiles = '/store/data/Run2012A/DoubleElectron/AOD/13Jul2012-v1/00000/00347915-EED9-E111-945A-848F69FD2817.root'
#options.inputFiles = '/store/user/andrey/Higgs_To_MuMuGamma_Dalitz_MH125_Mll_0to50_MadgraphHEFT_pythia6/AODSIM/39bf61f738ba3bdb8860f0848073cc88/aodsim_1_1_S9V.root'
#options.inputFiles = '/store/mc/Summer12_DR53X/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S10_START53_V7A-v1/0000/02CDCF05-BED2-E111-85F4-0030486740BA.root'
#options.inputFiles = '/store/data/Run2012D/SinglePhotonParked/AOD/22Jan2013-v1/30004/144D7268-4086-E211-9DC1-001E673984C1.root'
#options.inputFiles='/store/data/Run2012C/MuEG/AOD/22Jan2013-v1/20000/5A2C6379-8867-E211-BB9E-00266CFFC7E4.root'
#options.inputFiles = '/store/user/cmkuo/Dalitz_H_eeg_125_MG_v1/Dalitz_H_eeg_m125_RECO_v1/d459946fa1e058e24b305fca3ec661c6/STEP2_RAW2DIGI_L1Reco_RECO_PU_100_1_jL0.root'
#options.inputFiles = '/store/data/Run2012D/DoublePhoton/AOD/22Jan2013-v1/30002/FEA0AD2A-A284-E211-9B86-180373FF8D6A.root'
#options.inputFiles = '/store/data/Run2012C/DoubleMuParked/AOD/22Jan2013-v1/10000/00858723-296D-E211-A4B3-00259073E488.root'
#options.inputFiles = '/store/data/Run2012C/DoubleElectron/AOD/22Jan2013-v1/20000/00ABC56B-3668-E211-A0A5-003048678FDE.root'
#options.inputFiles = '/store/data/Run2011A/DoubleMu/AOD/21Jun2013-v1/10000/006AD75C-17DE-E211-B24E-003048678B0C.root'
options.inputFiles = '/store/mc/Fall13dr/JPsiToMuMu_Pt20to100-pythia6-gun/GEN-SIM-RECO/tsg_PU40bx50_POSTLS162_V2-v1/20000/06252816-B574-E311-B296-0025907DCA72.root'

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

# Basic MET shit
process.load("JetMETCorrections.Type1MET.correctionTermsPfMetType1Type2_cff")
process.load("JetMETCorrections.Type1MET.correctionTermsPfMetType0PFCandidate_cff")
process.load("JetMETCorrections.Type1MET.correctionTermsPfMetShiftXY_cff")
process.load("JetMETCorrections.Type1MET.correctedMet_cff")

if (isRealData):
  #process.GlobalTag.globaltag = 'FT_53_V21_AN3::All'
    process.GlobalTag.globaltag = 'FT_53_V21_AN4::All'
    process.corrPfMetType1.jetCorrLabel = cms.string("ak5PFL1FastL2L3Residual")
    process.corrPfMetShiftXY.parameter = process.pfMEtSysShiftCorrParameters_2012runABCDvsNvtx_data
else:
  #process.GlobalTag.globaltag = 'START53_V7N::All'
    process.GlobalTag.globaltag = 'START53_V23::All'
    process.corrPfMetType1.jetCorrLabel = cms.string("ak5PFL1FastL2L3")
    process.corrPfMetShiftXY.parameter = process.pfMEtSysShiftCorrParameters_2012runABCDvsNvtx_mc

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



if (isRealData):
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


process.kt6PFJetsIso = process.kt6PFJets.clone()
process.kt6PFJetsIso.doRhoFastjet = True
process.kt6PFJetsIso.Rho_EtaMax = cms.double(2.5)


# for jet pileup ID variables
from RecoJets.JetProducers.PileupJetIDParams_cfi import *

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
    fileName = cms.string('nuTuple.root')
    #fileName = cms.string('nuTuple_MuMu_Sync.root')
    #fileName = cms.string('nuTuple_EE_Sync.root')
    #fileName = cms.string('~/EOS/V09_05_8TeV/DoubleMu/nuTuple_Sync.root')
                                   )
# Electron Iso Stuff
# PF isolations for electronss
from CommonTools.ParticleFlow.Tools.pfIsolation import setupPFElectronIso
process.eleIsoSequence = setupPFElectronIso(process, 'gsfElectrons')



############################
### b-tag truth matching ###
############################

process.load("PhysicsTools.JetMCAlgos.CaloJetsMCFlavour_cfi")


# Flavour by reference
process.JetbyRef = cms.EDProducer("JetPartonMatcher",
                                  jets = cms.InputTag("ak5PFJetsL1FastL2L3"),
                                  coneSizeToAssociate = cms.double(0.3),
                                  partons = cms.InputTag("myPartons")
                                  )
# Flavour by value PhysDef
process.JetbyValPhys = cms.EDProducer("JetFlavourIdentifier",
                                      srcByReference = cms.InputTag("JetbyRef"),
                                      physicsDefinition = cms.bool(True),
                                      leptonInfo = cms.bool(True)
                                      )
# Flavour by value AlgoDef
process.JetbyValAlgo = cms.EDProducer("JetFlavourIdentifier",
                                      srcByReference = cms.InputTag("JetbyRef"),
                                      physicsDefinition = cms.bool(False),
                                      leptonInfo = cms.bool(True)
                                      )
process.JetFlavour = cms.Sequence(process.JetbyRef*process.JetbyValPhys*process.JetbyValAlgo)


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

# Photon CiC stuff
from NWU.ntupleProducer.hggPhotonIDCuts_cfi import *


### ntuple producer
process.ntupleProducer   = cms.EDAnalyzer('ntupleProducer',

  verboseTrigs         =    cms.untracked.bool(False),
  verboseMVAs          =    cms.untracked.bool(False),

  photonIsoCalcTag  =    cms.PSet(isolationSumsCalculator),
  jetPUIdAlgo       =    cms.PSet(full_5x),

  JetTag            =    cms.untracked.InputTag('ak5PFJetsL1FastL2L3'),
  JecTag            =    cms.string("AK5PF"),
  GenJetTag         =    cms.untracked.InputTag('ak5GenJets'),
  #METTag            =    cms.untracked.InputTag('pfMet'),
  METTag            =    cms.untracked.InputTag('pfMetT0pcT1Txy','','NTUPLE'),
  ElectronTag       =    cms.untracked.InputTag('gsfElectrons'),
  MuonTag           =    cms.untracked.InputTag('muons'),
  PhotonTag         =    cms.untracked.InputTag('photons'),
  PrimaryVtxTag     =    cms.untracked.InputTag('offlinePrimaryVertices'),
  rhoCorrTag        =    cms.untracked.InputTag('kt6PFJets',    'rho', 'RECO'),
  rho25CorrTag      =    cms.untracked.InputTag('kt6PFJetsIso', 'rho', 'NTUPLE'),
  rhoMuCorrTag      =    cms.untracked.InputTag('kt6PFJetsCentralNeutral', 'rho','RECO'),  # specifically for muon iso

  partFlowTag       =  cms.untracked.InputTag("particleFlow"), #,"Cleaned"),

  saveEleCrystals   =    cms.untracked.bool(False),
  savePhoCrystals   =    cms.untracked.bool(False),

  saveTriggerObj    =    cms.untracked.bool(False),

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

  isolation_cone_size_forSCremoval_pho = cms.untracked.double(0.3),
  isolation_cone_size_forSCremoval_el = cms.untracked.double(0.4),

  #for Ecal LazyTools and photon items
  ebReducedRecHitCollection = cms.InputTag("reducedEcalRecHitsEB"),
  eeReducedRecHitCollection = cms.InputTag("reducedEcalRecHitsEE"),
  esReducedRecHitCollection = cms.InputTag("reducedEcalRecHitsES"),
  hggPhotonIDConfiguration = hggPhotonIDCuts,

  #for proper electron pf iso calculation. EGamma group should be embarrassed that this is not in AOD
 # IsoDepElectron = cms.VInputTag(
 #   cms.InputTag('elPFIsoDepositChargedPFIso'),
 ##   cms.InputTag('elPFIsoDepositNeutralPFIso'),
 #   cms.InputTag('elPFIsoDepositGammaPFIso')),
 # IsoValElectronPF = cms.VInputTag(
 #   cms.InputTag('elPFIsoValueCharged03PFIdPFIso'),
 #   cms.InputTag('elPFIsoValueGamma03PFIdPFIso'),
 #   cms.InputTag('elPFIsoValueNeutral03PFIdPFIso')),

  #cms.InputTag('elPFIsoValueCharged04PFIdPFIso'),
  #  cms.InputTag('elPFIsoValueGamma04PFIdPFIso'),
   # cms.InputTag('elPFIsoValueNeutral04PFIdPFIso')),

  #Trigger stuff

  hltName           =    cms.untracked.string("HLT"),
  triggers          =    cms.untracked.vstring(
                                               "HLT_DoubleMu7_v",
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

                                               "HLT_Photon26_CaloId10_Iso50_Photon18_CaloId10_Iso50_Mass60_v",
                                               "HLT_Photon26_CaloId10_Iso50_Photon18_R9Id85_Mass60_v",
                                               "HLT_Photon26_R9Id85_Photon18_CaloId10_Iso50_Mass60_v",
                                               "HLT_Photon26_R9Id85_Photon18_R9Id85_Mass60_v",
                                               "HLT_Photon26_R9Id85_OR_CaloId10_Iso50_Photon18_R9Id85_OR_CaloId10_Iso50_Mass70_v",
                                               "HLT_Photon26_R9Id85_OR_CaloId10_Iso50_Photon18_R9Id85_OR_CaloId10_Iso50_Mass60_v",

                                               #some quarkonia triggers for h ->j/psi gamma analysis
                                               "HLT_Mu5_Track2_Jpsi_v",
                                               "HLT_Mu7_Track7_Jpsi_v",
                                               "HLT_Dimuon8_Jpsi_v",
                                               "HLT_Dimuon10_Jpsi_v",

                                               #"HLT_Photon30_R9Id90_CaloId_HE10_Iso40_EBOnly_Met25_HBHENoiseCleaned",
                                               #"HLT_Photon30_R9Id90_CaloId_HE10_Iso40_EBOnly",
                                               #"HLT_Photon30",
                                               "HLT_DiJet20_MJJ650_AllJets_DEta3p5_HT120_VBF",
                                               "HLT_DiJet30_MJJ700_AllJets_DEta3p5_VBF",
                                               "HLT_DiJet35_MJJ650_AllJets_DEta3p5_VBF",
                                               "HLT_DiJet35_MJJ700_AllJets_DEta3p5_VBF",
                                               "HLT_DiJet35_MJJ750_AllJets_DEta3p5_VBF"
                                               )
                                          )

process.eleIso = cms.Sequence(process.pfParticleSelectionSequence*process.eleIsoSequence)

process.preNtuple = cms.Sequence(
    process.goodOfflinePrimaryVertices
    * process.correctionTermsPfMetType1Type2
    * process.correctionTermsPfMetType0PFCandidate
    * process.correctionTermsPfMetShiftXY
    * process.pfMetT0pcT1Txy
    * process.kt6PFJetsIso
    * process.ak5PFJetsL1FastL2L3
    * process.ak5JetTracksAssociatorAtVertex
    * process.btagging
    * AllFilters

    * process.eleRegressionEnergy
    * process.calibratedElectrons
    * process.mvaTrigV0
    * process.mvaNonTrigV0

    * process.heepIdNoIso
    * process.heepIdNoIsoEles
    * process.modElectronIso
)

process.JetMC = cms.Sequence(process.myPartons * process.JetFlavour)

if isRealData:
    process.ntuplePath = cms.Path(process.preNtuple * process.eleIso * process.ntupleProducer)
else:
    process.ntuplePath = cms.Path(process.preNtuple * process.JetMC * process.eleIso * process.ntupleProducer)

#process.out = cms.OutputModule("PoolOutputModule",
#                               fileName = cms.untracked.string('/tmp/myTuple.root'),
#                               SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('ntuplePath')),
#                               outputCommands = cms.untracked.vstring('keep *')
#                               )
#process.outpath = cms.EndPath(process.out)

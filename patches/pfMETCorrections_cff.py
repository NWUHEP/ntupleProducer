import FWCore.ParameterSet.Config as cms

# load jet energy correction parameters
from JetMETCorrections.Configuration.JetCorrectionServices_cff import *

#--------------------------------------------------------------------------------
# select PFCandidates ("unclustered energy") not within jets
# for Type 2 MET correction
from CommonTools.ParticleFlow.TopProjectors.pfNoJet_cfi import pfNoJet
pfCandsNotInJet = pfNoJet.clone(
    topCollection = cms.InputTag('ak5PFJets'),
    bottomCollection = cms.InputTag('particleFlow')
)
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# produce Type 1 + 2 MET corrections for PFJets
pfJetMETcorr = cms.EDProducer("PFJetMETcorrInputProducer",
    src = cms.InputTag('ak5PFJets'),
    offsetCorrLabel = cms.string("ak5PFL1Fastjet"),
    jetCorrLabel = cms.string("ak5PFL1FastL2L3"), # NOTE: use "ak5PFL1FastL2L3" for MC / "ak5PFL1FastL2L3Residual" for Data
    jetCorrEtaMax = cms.double(9.9),
    type1JetPtThreshold = cms.double(10.0),
    skipEM = cms.bool(True),
    skipEMfractionThreshold = cms.double(0.90),
    skipMuons = cms.bool(True),
    skipMuonSelection = cms.string("isGlobalMuon | isStandAloneMuon")
)                                         
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# produce Type 2 MET corrections for selected PFCandidates
pfCandMETcorr = cms.EDProducer("PFCandMETcorrInputProducer",
    src = cms.InputTag('pfCandsNotInJet')                         
)   
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# produce Type 0 MET corrections for selected vertices
pfchsMETcorr = cms.EDProducer("PFchsMETcorrInputProducer",
    src = cms.InputTag('offlinePrimaryVertices'),
    goodVtxNdof = cms.uint32(4),
    goodVtxZ = cms.double(24)
)   
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# use MET corrections to produce Type 1 / Type 1 + 2 corrected PFMET objects
pfType1CorrectedMet = cms.EDProducer("CorrectedPFMETProducer",
    src = cms.InputTag('pfMet'),
    applyType0Corrections = cms.bool(False),
    srcCHSSums = cms.VInputTag(
        cms.InputTag('pfchsMETcorr', 'type0')
    ),
    type0Rsoft = cms.double(0.6),
    applyType1Corrections = cms.bool(True),
    srcType1Corrections = cms.VInputTag(
        cms.InputTag('pfJetMETcorr', 'type1')
    ),
    applyType2Corrections = cms.bool(False)
)   
#########################
#########################
pfType1CorrectedMetType0 = cms.EDProducer("CorrectedPFMETProducer",
    src = cms.InputTag('pfMet'),
    applyType0Corrections = cms.bool(False),
    srcCHSSums = cms.VInputTag(
        cms.InputTag('pfchsMETcorr', 'type0')
    ),
    type0Rsoft = cms.double(0.6),
    applyType1Corrections = cms.bool(True),
    srcType1Corrections = cms.VInputTag(
        cms.InputTag('pfMETcorrType0'),
        cms.InputTag('pfJetMETcorr', 'type1')
    ),
    applyType2Corrections = cms.bool(False)
)   

#########################
#########################

pfType1p2CorrectedMet = cms.EDProducer("CorrectedPFMETProducer",
    src = cms.InputTag('pfMet'),
    applyType0Corrections = cms.bool(False),
    srcCHSSums = cms.VInputTag(
        cms.InputTag('pfchsMETcorr', 'type0')
    ),
    type0Rsoft = cms.double(0.6),
    applyType1Corrections = cms.bool(True),
    srcType1Corrections = cms.VInputTag(
        cms.InputTag('pfJetMETcorr', 'type1')
    ),
    applyType2Corrections = cms.bool(True),
    srcUnclEnergySums = cms.VInputTag(
        cms.InputTag('pfJetMETcorr', 'type2'),
        cms.InputTag('pfJetMETcorr', 'offset'),
        cms.InputTag('pfCandMETcorr')                                    
    ),                              
    type2CorrFormula = cms.string("A"),
    type2CorrParameter = cms.PSet(
        A = cms.double(1.4)
    )
)   
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# define sequence to run all modules
producePFMETCorrections = cms.Sequence(
    pfCandsNotInJet
   * pfJetMETcorr
   * pfCandMETcorr
   * pfchsMETcorr
   * pfType1CorrectedMet
   * pfType1CorrectedMetType0
   * pfType1p2CorrectedMet
)
#--------------------------------------------------------------------------------

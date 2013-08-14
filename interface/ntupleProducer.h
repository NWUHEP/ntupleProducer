// system include files
#include <memory>
#include <string>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "DataFormats/GeometryVector/interface/LocalPoint.h"
#include "DataFormats/GeometryVector/interface/LocalVector.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/Math/interface/Vector3D.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/deltaR.h"

// Libraries for objects
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/JetReco/interface/JetID.h"

#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/METReco/interface/PFMETFwd.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/METCollection.h"
#include "DataFormats/METReco/interface/METFwd.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/Luminosity/interface/LumiSummary.h"
#include "DataFormats/Common/interface/ValueMap.h"

// Generator data formats
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/JetMatching/interface/JetFlavour.h"
#include "SimDataFormats/JetMatching/interface/JetFlavourMatching.h"
#include "SimDataFormats/JetMatching/interface/MatchedPartons.h"
#include "SimDataFormats/JetMatching/interface/JetMatchedPartons.h"

// PAT
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Tau.h"

//#include "RecoVertex/PrimaryVertexProducer/interface/VertexHigherPtSquared.h"

// JEC
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"

// Jet associators
#include "RecoJets/JetAssociationAlgorithms/interface/JetTracksAssociationDRCalo.h"
#include "RecoJets/JetAssociationAlgorithms/interface/JetTracksAssociationDRVertex.h"

// EGamma tools
#include "RecoEgamma/PhotonIdentification/interface/PhotonIsolationCalculator.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "EGamma/EGammaAnalysisTools/interface/ElectronEffectiveArea.h"
#include "EGamma/EGammaAnalysisTools/src/PFIsolationEstimator.cc"
//#include "EGamma/EGammaAnalysisTools/interface/EGammaMvaEleEstimator.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/METReco/interface/BeamHaloSummary.h"
#include "EGamma/EGammaAnalysisTools/interface/ElectronEnergyRegressionEvaluate.h"
#include "EgammaAnalysis/ElectronTools/interface/PatElectronEnergyCalibrator.h"

// Track MET Tools - Added by Rafael on May 28th
#include "RecoMET/METProducers/interface/TrackMETProducer.h"
#include "RecoMET/METProducers/interface/ParticleFlowForChargedMETProducer.h"

// ntuple storage classes
#include "TCPrimaryVtx.h"
#include "TCJet.h"
#include "TCMET.h"
#include "TCMuon.h"
#include "TCElectron.h"
#include "TCTau.h"
#include "TCPhoton.h"
#include "TCTriggerObject.h"
#include "TCGenJet.h"
#include "TCGenParticle.h"

// Need for HLT trigger info:
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"

//Root  stuff
#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TString.h"
#include "TObject.h"
#include "TObjArray.h"
#include "TClonesArray.h"
#include "TRefArray.h"
#include "TLorentzVector.h"
#include "TVector3.h"

using namespace edm;
using namespace std;
using namespace reco;

//
// class declaration
//

struct Filters {		//Filters 
  Bool_t isScraping;
  Bool_t isNoiseHcalHBHE;
  Bool_t isNoiseHcalLaser;
  Bool_t isNoiseEcalTP;
  Bool_t isNoiseEcalBE;
  Bool_t isCSCTightHalo;
  Bool_t isCSCLooseHalo;
  Bool_t isNoiseTracking;
  Bool_t isNoiseEEBadSc;
  Bool_t isNoisetrkPOG1;
  Bool_t isNoisetrkPOG2;
  Bool_t isNoisetrkPOG3;
};

/*
struct CrystalInfo{
  uint32_t rawId;
  int ieta;
  int iphi;
  int ix;
  int iy;
  double energy;
  double time;
  double timeErr;
  int recoFlag;
};
*/

class ntupleProducer : public EDAnalyzer {
 public:
  explicit ntupleProducer(const ParameterSet&);
  ~ntupleProducer();
  
 private:
  virtual void beginJob() ;
  virtual void beginRun(const Run&, const EventSetup&) ;
  virtual void analyze(const Event&, const EventSetup&);
  virtual void endLuminosityBlock(const LuminosityBlock&,const EventSetup&);
  virtual void endRun(const Run&, const EventSetup&);
  virtual void endJob() ;
  
  virtual bool  triggerDecision(Handle<TriggerResults>& hltR, int iTrigger);
  virtual float sumPtSquared(const Vertex& v);
  virtual bool  associateJetToVertex(reco::PFJet inJet, Handle<reco::VertexCollection> vtxCollection, TCJet *outJet);   
  virtual void  electronMVA(const reco::GsfElectron* iElectron, TCElectron* eleCon, const Event& iEvent,const EventSetup& iSetup, const reco::PFCandidateCollection& PFCandidates, float Rho);
  virtual bool  isFilteredOutScraping(const Event& iEvent, const EventSetup& iSetup, int numtrack=10, double thresh=0.25);
  virtual float MatchBTagsToJets(const reco::JetTagCollection, const reco::PFJet);
  void analyzeTrigger(Handle<TriggerResults> &hltR, Handle<trigger::TriggerEvent> &hltE, const std::string& triggerName, int* trigCount);                   
  // ----------member data ---------------------------
  
  struct JetCompare :
    public std::binary_function<reco::Jet, reco::Jet, bool> {
    inline bool operator () (const reco::Jet &j1,
			     const reco::Jet &j2) const
    { return (j1.p4().Pt() > j2.p4().Pt()); }
  };
  
    static bool EnergySortCriterium( const TCPhoton::CrystalInfo p1,const TCPhoton::CrystalInfo p2 ){
      return p1.energy > p2.energy;
    };
    
  typedef std::map<reco::Jet, unsigned int, JetCompare> flavourMap;
  typedef reco::JetTagCollection::const_iterator tag_iter;
  
  //Standard event info
  ULong64_t   eventNumber;
  UInt_t      runNumber, lumiSection, bunchCross, nEvents;
  float ptHat, qScale, evtWeight;
  float deliveredLumi, recordedLumi, lumiDeadTime;
  float rhoFactor, rho25Factor, rhoMuFactor,;
  vector<string>  savedTriggerNames;
  
  Service<TFileService> fs;
  TTree* eventTree;
  TTree* jobTree;

  InputTag jetTag_;
  InputTag metTag_;
  InputTag genJetTag_;
  InputTag muonTag_;
  InputTag electronTag_;
  InputTag photonTag_;
  InputTag primaryVtxTag_;
  InputTag triggerResultsTag_;
  InputTag rhoCorrTag_, rho25CorrTag_, rhoMuCorrTag_;
  InputTag hcalHBHEFilterTag_;
  InputTag hcalLaserFilterTag_;
  InputTag trackingFailureTag_;
  InputTag eeBadScFilterTag_;
  InputTag trkPOGFiltersTag1_;
  InputTag trkPOGFiltersTag2_;
  InputTag trkPOGFiltersTag3_;
  InputTag partFlowTag_;
  ParameterSet photonIsoCalcTag_;
  InputTag triggerEventTag_;

  bool saveJets_;
  bool saveElectrons_;
  bool saveMuons_;
  bool savePhotons_;
  bool saveMET_;
  bool saveGenJets_;
  bool saveGenParticles_;
  bool isRealData;
  bool verboseTrigs;
  bool verboseMVAs;
  
  //Physics object containers
  TClonesArray* recoJets;
  TClonesArray* recoMuons;
  TClonesArray* recoElectrons;
  TClonesArray* recoPhotons;
  TClonesArray* triggerObjects;
  TClonesArray* genJets;
  TClonesArray* genParticles;
  TCMET*        recoMET;
  TCMET*	    track_MET; 
  TCMET*	    T0MET; 
  TCMET*	    T2MET; 
  TCMET*        recoMET_corr;
  TCMET*        mva_MET;

  //Vertex info
  TClonesArray* primaryVtx;
  TVector3*     beamSpot;
  unsigned      nPUVertices;
  float         nPUVerticesTrue;
  
  //Triggers
  HLTConfigProvider hltConfig_;
  string            hlTriggerResults_, hltProcess_, triggerName_;
  TriggerNames      triggerNames;
  vector<string>    hlNames;
  vector<string>    triggerPaths_;
  ULong64_t         triggerStatus;
  unsigned          hltPrescale[64];
  
  // Technical filters
  Filters myNoiseFilters;

  //Isolator
  PFIsolationEstimator phoIsolator;
  PFIsolationEstimator eleIsolator;

  // Electron Regression
  ElectronEnergyRegressionEvaluate* myEleReg;
  
};

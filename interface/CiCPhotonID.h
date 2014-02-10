//
// Original Author:  Matteosan SANI
// Created:  Tue Oct  18 10:14:43 CET 2011
// 
//
#ifndef CICPHOTONID_H
#define CICPHOTONID_H

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "TLorentzVector.h"

class CiCPhotonID {
 public:
  
  enum CiCPhotonIDLevel {NoCut = -1, Loose, Medium, Tight, SuperTight, HyperTight1, HyperTight2, HyperTight3, HyperTight4, 
       HyperTight5, HyperTight6, HyperTight7, Preselection};
  
  explicit CiCPhotonID(const edm::ParameterSet&);
  ~CiCPhotonID() {};

  TLorentzVector get_pho_p4(reco::PhotonRef, Int_t);
  int PhotonIDCategory(reco::PhotonRef, int);
  Float_t DeltaRToTrackHgg(reco::PhotonRef, Int_t);
  bool PhotonIDPF(int, reco::PhotonRef, Int_t, CiCPhotonIDLevel);
  bool PhotonID(int, reco::PhotonRef, int, CiCPhotonIDLevel);
  Float_t SumTrackPtInConeHgg(reco::PhotonRef, Int_t, Float_t, Float_t, Float_t, Float_t, Float_t, Float_t);
  Float_t WorstSumTrackPtInConeHgg(reco::PhotonRef, Float_t, Float_t, Float_t, Float_t, Float_t, Float_t);

  int photonCutLevel4cat(reco::PhotonRef, Int_t);
  int photonCutLevel6cat(reco::PhotonRef, Int_t);
  int photonCutLevel6catPF(reco::PhotonRef, Int_t);

  std::vector<float> pfTkIsoWithVertex(reco::PhotonRef, float, float, float, float, float, float, reco::PFCandidate::ParticleType);
  float pfEcalIso(reco::PhotonRef, float, float, float, float, float, float, float, reco::PFCandidate::ParticleType);
  float pfHcalIso(reco::PhotonRef, float, float, reco::PFCandidate::ParticleType);
  
  void setPhotonIDThresholds(const edm::ParameterSet&);
  void configure(edm::Handle<reco::VertexCollection>, 
     edm::Handle<reco::TrackCollection>,
     edm::Handle<reco::GsfElectronCollection>, 
     edm::Handle<reco::PFCandidateCollection>,
     double);
  
 private:
  edm::Handle<reco::VertexCollection> vtxHandle;
  edm::Handle<reco::TrackCollection> tkHandle;
  edm::Handle<reco::GsfElectronCollection> elHandle;
  edm::Handle<reco::PFCandidateCollection> pfHandle;
  edm::Handle<double> rhoHandle;

  double rho; 

  // PhotonID thresholds
  std::vector<double> phoIDcuts[12][7];
  std::vector<double> phoIDcuts6cat[12][7];
  std::vector<double> phoIDcuts6catPF[12][7];
};

#endif

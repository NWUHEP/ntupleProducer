#define MSDEBUG 0
#include "../interface/CiCPhotonID.h"

#include "DataFormats/Math/interface/deltaR.h"

CiCPhotonID::CiCPhotonID(const edm::ParameterSet& config) {
  setPhotonIDThresholds(config);
}

void CiCPhotonID::configure(edm::Handle<reco::VertexCollection> vtxs, edm::Handle<reco::TrackCollection> tks,
          edm::Handle<reco::GsfElectronCollection> els, edm::Handle<reco::PFCandidateCollection> pfs,
          double myRho) {

  vtxHandle = vtxs;
  tkHandle = tks;
  elHandle = els;
  pfHandle = pfs;
  rho = myRho;
}

TLorentzVector CiCPhotonID::get_pho_p4(reco::PhotonRef photon, Int_t ivtx) {
  
  reco::VertexRef v(vtxHandle, ivtx);
  
  TVector3 vtx = TVector3(v->x(), v->y(), v->z());
  TVector3 pho = TVector3(photon->caloPosition().x(), photon->caloPosition().y(), photon->caloPosition().z());
  TVector3 direction = pho - vtx;
  TVector3 p = direction.Unit() * photon->energy();
  TLorentzVector p4(p.x(), p.y(), p.z(), photon->energy());

  return p4;
}

int CiCPhotonID::PhotonIDCategory(reco::PhotonRef photon, int nCategories) {
  
  float thisCat = -1; 

  float eta = fabs(photon->caloPosition().Eta());
  float r9  = photon->r9();
  
  if (nCategories == 4) {
    bool etaCat = eta>1.479;
    bool r9Cat  = r9<0.94;
    
    thisCat = 2*(etaCat) + r9Cat;
  } else if (nCategories == 6) {
    Int_t r9Cat = (Int_t)(r9<0.94) + (Int_t)(r9<0.9); // 0, 1, or 2 (high r9 --> low r9)
    bool etaCat = eta>1.479;

    thisCat = 3*(etaCat) + r9Cat;
  } else {
    std::cout << "Wrong number of categories." << std::endl;
  }
    
  return thisCat;
}

Float_t CiCPhotonID::DeltaRToTrackHgg(reco::PhotonRef photon, Int_t maxlosthits) {

  if (MSDEBUG)
    std::cout << "DeltaRToTrackHgg BEGIN" << std::endl;

  int elind = -1;
  float eldr = 99.;

  for(unsigned int iel = 0; iel<elHandle->size(); ++iel) {

    reco::GsfElectronRef el(elHandle, iel);
    
    if(el->gsfTrack()->trackerExpectedHitsInner().numberOfHits() > maxlosthits)
      continue;

    if(el->superCluster() == photon->superCluster()) {
      elind = iel;    
      break;
    }
  }

  if(elind >= 0) {
    reco::GsfElectronRef el(elHandle, elind);
    eldr = sqrt(pow(el->deltaEtaSuperClusterTrackAtVtx(), 2) + pow(el->deltaPhiSuperClusterTrackAtVtx(), 2));
  }

  return eldr;
}

bool CiCPhotonID::PhotonIDPF(int nCategories, reco::PhotonRef photon, Int_t ivtx, CiCPhotonID::CiCPhotonIDLevel IDlevel) {

  //IDlevel 0 nocut, 1 Loose, 2 Medium...
  if (IDlevel == -1)
    return true;
  
  bool passes = true;
  float vars[10];
  
  int thisCat = PhotonIDCategory(photon, nCategories);
  
  TLorentzVector phop4 = get_pho_p4(photon, ivtx);
  
  float val_pfiso_photon04 = pfEcalIso(photon, 0.4, 0.045, 0.070, 0.015, 0.0, 0.08, 0.1, reco::PFCandidate::gamma); 
  float val_pfiso_photon03 = pfEcalIso(photon, 0.3, 0.045, 0.070, 0.015, 0.0, 0.08, 0.1, reco::PFCandidate::gamma);

  std::vector<float> vtxIsolations03 = pfTkIsoWithVertex(photon, 0.3, 0.02, 0.02, 0.0, 0.2, 0.1, reco::PFCandidate::h);
  std::vector<float> vtxIsolations04 = pfTkIsoWithVertex(photon, 0.4, 0.02, 0.02, 0.0, 0.2, 0.1, reco::PFCandidate::h);

  float val_pfiso_charged03 = vtxIsolations03[ivtx];
  float val_pfiso_charged_badvtx_04 = -99;
  //int badind = -1;
  for(unsigned int iv=0; iv<vtxHandle->size(); iv++) {
    if(vtxIsolations04[iv] > val_pfiso_charged_badvtx_04) {
      //badind = iv;
      val_pfiso_charged_badvtx_04 = vtxIsolations04[iv];
    }
  }

  float rhofacpf[6]    = {0.075, 0.082, 0.143, 0.050, 0.091, 0.106};
  float rhofacbadpf[6] = {0.141, 0.149, 0.208, 0.135, 0.162, 0.165};
  float rhofac         = rhofacpf[thisCat];
  float rhofacbad      = rhofacbadpf[thisCat];
  float isosumconst    = 2.8;
  float isosumconstbad = 4.8;

  float val_isosum     = val_pfiso_charged03 + val_pfiso_photon03;
  float val_isosumbad  = val_pfiso_charged_badvtx_04 + val_pfiso_photon04;
  float val_isosumoet    = (val_isosum + isosumconst - rho*rhofac)*50./phop4.Et();
  float val_isosumoetbad = (val_isosumbad + isosumconstbad - rho*rhofacbad)*50./phop4.Et();
  float val_trkisooet    = (val_pfiso_charged03)*50./phop4.Et();

  vars[0] = val_isosumoet;
  vars[1] = val_isosumoetbad;
  vars[2] = val_trkisooet;
  vars[3] = photon->sigmaIetaIeta();
  vars[4] = photon->hadronicOverEm();
  vars[5] = photon->r9();
  vars[6] = DeltaRToTrackHgg(photon,  99);
  //vars[7] = (float)photon->hasPixelSeed();
  
  for(int var=0; var<7; var++){

    if (nCategories == 6) {
      if(var == 5 || var == 6){   //cuts from below
  if(vars[var] < phoIDcuts6catPF[IDlevel][var][thisCat]) {
    passes = false;
    break;
  }
      } else {                    //cuts from above
  if (vars[var] > phoIDcuts6catPF[IDlevel][var][thisCat]) {
    passes = false;
    break;
  }
      }
    }
  } 

  //std::cout << passes << std::endl;
  return passes;
}

bool CiCPhotonID::PhotonID(int nCategories, reco::PhotonRef photon, int iVtx, CiCPhotonID::CiCPhotonIDLevel IDlevel) {

 //IDlevel 0 nocut, 1 Loose, 2 Medium...
  if ((int)IDlevel == -1)
    return true;
  

  bool passes = true;
  float vars[10];
  float rhofacbad=0.52, rhofac=0.17;
  
  int thisCat = PhotonIDCategory(photon, nCategories);
  TLorentzVector phop4 = get_pho_p4(photon, iVtx);

  float pho_tkiso_goodvtx = SumTrackPtInConeHgg(photon, iVtx, 0, 0.30, 0.02, 0.0, 1.0, 0.1);
  float pho_tkiso_badvtx = WorstSumTrackPtInConeHgg(photon, 0, 0.40, 0.02, 0.0, 1.0, 0.1);
  float isosumoet = (pho_tkiso_goodvtx + photon->ecalRecHitSumEtConeDR03() + 
         photon->hcalTowerSumEtConeDR04() - rho*rhofac)*50./phop4.Et();
  float isosumoetbad = (pho_tkiso_badvtx + photon->ecalRecHitSumEtConeDR04() +
      photon->hcalTowerSumEtConeDR04() - rho*rhofacbad)*50./phop4.Et();
  float trkisooet = (pho_tkiso_goodvtx)*50./phop4.Et();

  vars[0] = isosumoet;
  vars[1] = isosumoetbad;
  vars[2] = trkisooet;
  vars[3] = photon->sigmaIetaIeta();
  vars[4] = photon->hadronicOverEm();
  vars[5] = photon->r9();
  vars[6] = DeltaRToTrackHgg(photon, 99);
  //vars[7] = (float)photon->hasPixelSeed();
  
  for(int var=0; var<7; var++){
    if (nCategories == 4) {
      if(var == 5 || var == 6){   //cuts from below
  if(vars[var] < phoIDcuts[IDlevel][var][thisCat]) {
    passes =false;
    break;
  }
      } else {                    //cuts from above
  if (vars[var] > phoIDcuts[IDlevel][var][thisCat]) {
    passes =false;
    break;
  }
      }
    }

    if (nCategories == 6) {
      if(var == 5 || var == 6){   //cuts from below
  if(vars[var] < phoIDcuts6cat[IDlevel][var][thisCat]) {
    passes =false;
    break;
  }
      } else {                    //cuts from above
  if (vars[var] > phoIDcuts6cat[IDlevel][var][thisCat]) {
    passes =false;
    break;
  }
      }
    }
  } 

  return passes;
}

Float_t CiCPhotonID::SumTrackPtInConeHgg(reco::PhotonRef photon, Int_t iVtx, Float_t PtMin, Float_t OuterConeRadius, Float_t InnerConeRadius, Float_t EtaStripHalfWidth, Float_t dzmax, Float_t dxymax) {

  TLorentzVector phop4 = get_pho_p4(photon, iVtx);

  reco::VertexRef v(vtxHandle, iVtx);
  TVector3 vtxpos = TVector3(v->x(), v->y(), v->z());
  
  float SumTrackPt=0;
  for(unsigned int itk=0; itk<tkHandle->size(); itk++) {
    
    reco::TrackRef tk(tkHandle, itk);
    
    if(tk->pt() < PtMin)
      continue;
    
    double deltaz = fabs(v->z() - tk->vz());
    if(deltaz > dzmax)
      continue;

    double dxy = (-(tk->vx() - v->x())*tk->py() + (tk->vy() - v->y())*tk->px()) / tk->pt();
    if(fabs(dxy) > dxymax)
      continue;

    double tk_eta = tk->eta();
    double tk_phi = tk->phi();
    double deta = fabs(phop4.Eta() - tk_eta);
    double dphi = fabs(phop4.Phi() - tk_phi);
    if(dphi > TMath::Pi())
      dphi = TMath::TwoPi() - dphi;

    double deltaR = sqrt(deta*deta + dphi*dphi);
    if(deltaR < OuterConeRadius && deltaR >= InnerConeRadius && deta >= EtaStripHalfWidth)
      SumTrackPt += tk->pt();
  }

  return SumTrackPt;
}

Float_t CiCPhotonID::WorstSumTrackPtInConeHgg(reco::PhotonRef photon, Float_t PtMin, Float_t OuterConeRadius, Float_t InnerConeRadius, Float_t EtaStripHalfWidth, Float_t dzmax, Float_t dxymax) {

  //Int_t worstvtxind = -1;
  Float_t maxisosum = -100;
  for(unsigned int ivtx=0; ivtx<vtxHandle->size(); ++ivtx) {
    Float_t thisvtxisosum = SumTrackPtInConeHgg(photon, ivtx, PtMin, OuterConeRadius, InnerConeRadius, EtaStripHalfWidth, dzmax, dxymax);
    if(thisvtxisosum > maxisosum) {
      maxisosum = thisvtxisosum;
      //worstvtxind = ivtx;
    }
  }

  return maxisosum;
}

std::vector<float> CiCPhotonID::pfTkIsoWithVertex(reco::PhotonRef localPho, float dRmax, float dRvetoBarrel, float dRvetoEndcap, 
              float ptMin, float dzMax, float dxyMax,
              reco::PFCandidate::ParticleType pfToUse) {
  
  float dRveto;
  if (localPho->isEB())
    dRveto = dRvetoBarrel;
  else
    dRveto = dRvetoEndcap;

  std::vector<float> result;
  const reco::PFCandidateCollection* forIsolation = pfHandle.product();

  //Calculate isolation sum separately for each vertex
  for(unsigned int ivtx=0; ivtx<vtxHandle->size(); ++ivtx) {
    
    // Shift the photon according to the vertex
    reco::VertexRef vtx(vtxHandle, ivtx);
    math::XYZVector photon_directionWrtVtx(localPho->superCluster()->x() - vtx->x(),
             localPho->superCluster()->y() - vtx->y(),
             localPho->superCluster()->z() - vtx->z());
    
    float sum = 0;
    // Loop over the PFCandidates
    for(unsigned i=0; i<forIsolation->size(); i++) {
      
      const reco::PFCandidate& pfc = (*forIsolation)[i];

      //require that PFCandidate is a charged hadron
      if (pfc.particleId() == pfToUse) {
  if (pfc.pt() < ptMin)
    continue;

  float dz = fabs(pfc.trackRef()->dz(vtx->position()));
  if (dz > dzMax) continue;

  float dxy = fabs(pfc.trackRef()->dxy(vtx->position()));
  if(fabs(dxy) > dxyMax) continue;

  float dR = deltaR(photon_directionWrtVtx.Eta(), photon_directionWrtVtx.Phi(), pfc.momentum().Eta(), pfc.momentum().Phi());
  if(dR > dRmax || dR < dRveto) continue;

  sum += pfc.pt();
      }
    }

    result.push_back(sum);
  }
  
  return result;
}

float CiCPhotonID::pfEcalIso(reco::PhotonRef localPho, float dRmax, float dRVetoBarrel, float dRVetoEndcap, float etaStripBarrel, float etaStripEndcap, float energyBarrel, float energyEndcap, reco::PFCandidate::ParticleType pfToUse) {
  
  float dRVeto;
  float etaStrip;

  if (localPho->isEB()) {
    dRVeto = dRVetoBarrel;
    etaStrip = etaStripBarrel;
  } else {
    dRVeto = dRVetoEndcap;
    etaStrip = etaStripEndcap;
  }
      
  const reco::PFCandidateCollection* forIsolation = pfHandle.product();

  float sum = 0;
  for(unsigned i=0; i<forIsolation->size(); i++) {
    
    const reco::PFCandidate& pfc = (*forIsolation)[i];
    
    if (pfc.particleId() ==  pfToUse) {
      
      // Do not include the PFCandidate associated by SC Ref to the reco::Photon
      if(pfc.superClusterRef().isNonnull() && localPho->superCluster().isNonnull()) {
  if (pfc.superClusterRef() == localPho->superCluster()) 
    continue;
      }
      
      if (localPho->isEB()) {
  if (fabs(pfc.pt()) < energyBarrel)
    continue;
      } else {
  if (fabs(pfc.energy()) < energyEndcap)
    continue;
      }
      
      // Shift the photon direction vector according to the PF vertex
      math::XYZPoint pfvtx = pfc.vertex();
      math::XYZVector photon_directionWrtVtx(localPho->superCluster()->x() - pfvtx.x(),
               localPho->superCluster()->y() - pfvtx.y(),
               localPho->superCluster()->z() - pfvtx.z());
    
      float dEta = fabs(photon_directionWrtVtx.Eta() - pfc.momentum().Eta());
      float dR = deltaR(photon_directionWrtVtx.Eta(), photon_directionWrtVtx.Phi(), pfc.momentum().Eta(), pfc.momentum().Phi());
      
      if (dEta < etaStrip)
  continue;
      
      if(dR > dRmax || dR < dRVeto)
  continue;
      
      sum += pfc.pt();
    }
  }
  
  return sum;
}

float CiCPhotonID::pfHcalIso(reco::PhotonRef localPho, float dRmax, float dRveto, reco::PFCandidate::ParticleType pfToUse) {
  
  return pfEcalIso(localPho, dRmax, dRveto, dRveto, 0.0, 0.0, 0.0, 0.0, pfToUse);
}

void CiCPhotonID::setPhotonIDThresholds(const edm::ParameterSet& iConfig) {

  const edm::ParameterSet gammaIDCuts = iConfig.getParameter<edm::ParameterSet>("hggPhotonIDConfiguration") ;
  char a[100];
  
  for (int lev = 0; lev < 12; ++lev) {
    sprintf(a, "cutsubleadisosumoet6c%d", lev);
    phoIDcuts6cat[lev][0]     = gammaIDCuts.getParameter<std::vector<double> >(a);
    sprintf(a, "cutsubleadisosumoetbad6c%d", lev);
    phoIDcuts6cat[lev][1]  = gammaIDCuts.getParameter<std::vector<double> >(a);
    sprintf(a, "cutsubleadtrkisooetom6c%d", lev);
    phoIDcuts6cat[lev][2]   = gammaIDCuts.getParameter<std::vector<double> >(a);
    sprintf(a, "cutsubleadsieie6c%d", lev);
    phoIDcuts6cat[lev][3]        = gammaIDCuts.getParameter<std::vector<double> >(a);
    sprintf(a, "cutsubleadhovere6c%d", lev);
    phoIDcuts6cat[lev][4]       = gammaIDCuts.getParameter<std::vector<double> >(a);
    sprintf(a, "cutsubleadr96c%d", lev);
    phoIDcuts6cat[lev][5]          = gammaIDCuts.getParameter<std::vector<double> >(a);
    sprintf(a, "cutsublead_drtotk_25_996c%d", lev);
    phoIDcuts6cat[lev][6] = gammaIDCuts.getParameter<std::vector<double> >(a);
    
    sprintf(a, "cutsubleadisosumoet%d", lev);
    phoIDcuts[lev][0]     = gammaIDCuts.getParameter<std::vector<double> >(a);
    sprintf(a, "cutsubleadisosumoetbad%d", lev);
    phoIDcuts[lev][1]  = gammaIDCuts.getParameter<std::vector<double> >(a);
    sprintf(a, "cutsubleadtrkisooetom%d", lev);
    phoIDcuts[lev][2]   = gammaIDCuts.getParameter<std::vector<double> >(a);
    sprintf(a, "cutsubleadsieie%d", lev);
    phoIDcuts[lev][3]        = gammaIDCuts.getParameter<std::vector<double> >(a);
    sprintf(a, "cutsubleadhovere%d", lev);
    phoIDcuts[lev][4]        = gammaIDCuts.getParameter<std::vector<double> >(a);
    sprintf(a, "cutsubleadr9%d", lev);
    phoIDcuts[lev][5]           = gammaIDCuts.getParameter<std::vector<double> >(a);
    sprintf(a, "cutsublead_drtotk_25_99%d", lev);
    phoIDcuts[lev][6] = gammaIDCuts.getParameter<std::vector<double> >(a);

    sprintf(a, "cutsubleadpfisosumoet%d", lev);
    phoIDcuts6catPF[lev][0]     = gammaIDCuts.getParameter<std::vector<double> >(a);
    sprintf(a, "cutsubleadpfisosumoetbad%d", lev);
    phoIDcuts6catPF[lev][1]  = gammaIDCuts.getParameter<std::vector<double> >(a);
    sprintf(a, "cutsubleadpfchisooet%d", lev);
    phoIDcuts6catPF[lev][2]   = gammaIDCuts.getParameter<std::vector<double> >(a);
    sprintf(a, "cutsubleadpfsieie%d", lev);
    phoIDcuts6catPF[lev][3]        = gammaIDCuts.getParameter<std::vector<double> >(a);
    sprintf(a, "cutsubleadpfhovere%d", lev);
    phoIDcuts6catPF[lev][4]       = gammaIDCuts.getParameter<std::vector<double> >(a);
    sprintf(a, "cutsubleadpfr9%d", lev);
    phoIDcuts6catPF[lev][5]          = gammaIDCuts.getParameter<std::vector<double> >(a);
    sprintf(a, "cutsubleadpf_drtotk_25_99%d", lev);
    phoIDcuts6catPF[lev][6] = gammaIDCuts.getParameter<std::vector<double> >(a);

  }
}

int CiCPhotonID::photonCutLevel4cat(reco::PhotonRef photon, Int_t iVtx) {
  
  int lev = 0;

  for (Int_t i=0; i<8; i++)
    if (PhotonID(4, photon, iVtx, CiCPhotonID::CiCPhotonIDLevel(i)))
      lev++;

  return lev;
}


int CiCPhotonID::photonCutLevel6cat(reco::PhotonRef photon, Int_t iVtx) {
  
  int lev = 0;

  for (Int_t i=0; i<8; i++)
    if (PhotonID(6, photon, iVtx, CiCPhotonID::CiCPhotonIDLevel(i)))
      lev++;

  return lev;
}

int CiCPhotonID::photonCutLevel6catPF(reco::PhotonRef photon, Int_t iVtx) {
  
  int lev = 0;

  for (Int_t i=0; i<11; i++)
    if (PhotonIDPF(6, photon, iVtx, CiCPhotonID::CiCPhotonIDLevel(i)))
      lev++;

  return lev;
}



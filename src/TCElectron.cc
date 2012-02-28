#include "Higgs/ntupleProducer/interface/TCElectron.h"
#include <iostream>

TCElectron::TCElectron() {
}

TCElectron::~TCElectron() {
}



// "get" methods -------------------------------------

TLorentzVector TCElectron::P4() const {
   return _p4;
}

TVector2 TCElectron::P2() const {
  TVector2 v2(_p4.Px(), _p4.Py());
  return v2;
}

float TCElectron::Pt() const {
   return _p4.Pt();
}

TVector3 TCElectron::Vtx() const {
   return _vtx;
}

float TCElectron::PtError() const {
  return _ptError;
}

float TCElectron::Et() const {
  return _p4.Et();
}

float TCElectron::Eta() const {
 return _p4.Eta();
}

float TCElectron::Phi() const {
  return _p4.Phi();
}

int TCElectron::Charge() const {
   return _charge;
}

int TCElectron::NumberOfValidPixelHits() const {
  return _numberOfValidPixelHits;
}

int TCElectron::NumberOfValidTrackerHits() const {
  return _numberOfValidTrackerHits;
}

int TCElectron::NumberOfLostPixelHits() const {
  return _numberOfLostPixelHits;
}

int TCElectron::NumberOfLostTrackerHits() const {
  return _numberOfLostTrackerHits;
}

float TCElectron::NormalizedChi2() const {
  return _normalizedChi2;
}


float TCElectron::EmIso() const {
   return _emIso03;
}
float TCElectron::HadIso() const {
   return _hadIso03;
}
float TCElectron::TrkIso() const {
   return _trkIso03;
}

float TCElectron::EmIso03() const {
   return _emIso03;
}
float TCElectron::HadIso03() const {
   return _hadIso03;
}
float TCElectron::TrkIso03() const {
   return _trkIso03;
}

float TCElectron::EmIso04() const {
   return _emIso04;
}
float TCElectron::HadIso04() const {
   return _hadIso04;
}
float TCElectron::TrkIso04() const {
   return _trkIso04;
}

float TCElectron::PfRelIso(float coneSize) const {
  float relIso = 0;
  if (fabs(coneSize - 0.3) < 0.01)
    relIso = (_pfIso_Pt03 + _pfIso_Gamma03 + _pfIso_Neutral03) / _p4.Pt();
  if (fabs(coneSize - 0.4) < 0.01)
    relIso = (_pfIso_Pt04 + _pfIso_Gamma04 + _pfIso_Neutral04) / _p4.Pt();
  if (fabs(coneSize - 0.5) < 0.01)
    relIso = (_pfIso_Pt05 + _pfIso_Gamma05 + _pfIso_Neutral05) / _p4.Pt();
  return relIso;
}

float TCElectron::PfSumPt(float coneSize) const {
  float sumPt = 0;
  if (fabs(coneSize - 0.3) < 0.01) sumPt = _pfIso_Pt03;
  if (fabs(coneSize - 0.4) < 0.01) sumPt = _pfIso_Pt04;
  if (fabs(coneSize - 0.5) < 0.01) sumPt = _pfIso_Pt05;
  return sumPt;
}

float TCElectron::PfENeutral(float coneSize) const {
  float neutral = 0;
  if (fabs(coneSize - 0.3) < 0.01) neutral = _pfIso_Neutral03;
  if (fabs(coneSize - 0.4) < 0.01) neutral = _pfIso_Neutral04;
  if (fabs(coneSize - 0.5) < 0.01) neutral = _pfIso_Neutral05;
  return neutral;
}

float TCElectron::PfEGamma(float coneSize) const {
  float gamma = 0;
  if (fabs(coneSize - 0.3) < 0.01) gamma = _pfIso_Gamma03;
  if (fabs(coneSize - 0.4) < 0.01) gamma = _pfIso_Gamma04;
  if (fabs(coneSize - 0.5) < 0.01) gamma = _pfIso_Gamma05;
  return gamma;
}

bool TCElectron::IsEB() const {
  return _isEB;
}

bool TCElectron::IsEE() const {
  return _isEE;
}

bool TCElectron::IsInGap() const {
  return _isInGap;
}


float TCElectron::HadOverEm() const {
  return _hadOverEm;
}
float TCElectron::DphiSuperCluster() const {
  return _dPhiSuperCluster;
}
float TCElectron::DetaSuperCluster() const {
  return _dEtaSuperCluster;
}
float TCElectron::SigmaIetaIeta() const {
  return _sigmaIetaIeta;
}

int TCElectron::ConversionFlag() const {
  return _convFlag;
}

float TCElectron::ConversionDist() const {
  return _convDist;
}

float TCElectron::ConversionDcot() const {
  return _convDcot;
}

float TCElectron::ConversionRad() const {
  return _convRad;
}

int TCElectron::CutLevel(int lvl) const{
  if(lvl==95){
    return _cut95;
  }else if(lvl==90) {
    return _cut90;
  }else if(lvl==85) {
    return _cut85;
  }else if(lvl==80) {
    return _cut80;
  }else if(lvl==70) {
    return _cut70;
  }else if(lvl==60) {
    return _cut60;
  }else{
    return -99;
  }
}

bool TCElectron::PassID(int lvl) const { 
  unsigned c = CutLevel(lvl);

  if (c & 0x01) return true;
  else return false;

  //if( (c == 1) || (c == 3) || (c == 5) || (c == 7)) return true;
  //return false;
}   

bool TCElectron::PassConversion(int lvl) const {
  int c = CutLevel(lvl);
  if(c < 0) return false;
  if( (c == 4) || (c == 5) || (c == 6) || (c == 7)) return true;
  return false;
}

bool TCElectron::PassIsolation(int lvl) const {
  int c = CutLevel(lvl);
  if(c < 0) return false;
  if( (c == 2) || (c == 3) || (c == 6) || (c == 7)) return true;
  return false;
}



float TCElectron::Dxy(TVector3 *primVtx) const {
  //Calculating track dxy parameter wrt primary vertex
  //d0 = - dxy
  float vx = _vtx.X(), vy = _vtx.Y();
  float px = _p4.Px(), py = _p4.Py(), pt = _p4.Pt();
  float pvx = primVtx->X(), pvy = primVtx->Y();
  float ret =  (-(vx-pvx)*py + (vy-pvy)*px)/pt;
  return ret;
}

float TCElectron::Dz(TVector3 *primVtx) const {
  //Calculating track dz parameter wrt primary vertex
  float vx = _vtx.X(), vy = _vtx.Y(), vz = _vtx.Z();
  float px = _p4.Px(), py = _p4.Py();
  float pz = _p4.Pz(), pt = _p4.Pt();
  float pvx = primVtx->X(), pvy = primVtx->Y(), pvz = primVtx->Z();
  float ret =  (vz-pvz)-((vx-pvx)*px +(vy-pvy)*py)/pt*(pz/pt);
  return ret;
}

//------------------------------------------------
// "set" methods ---------------------------------------------
//------------------------------------------------------------------------

void TCElectron::SetP4(TLorentzVector p4) {
  _p4 = p4;
}

void TCElectron::SetP4(float px, float py, float pz, float e) {
  TLorentzVector p4(px, py, pz, e);
  _p4 = p4;
}

void TCElectron::SetVtx(float vx, float vy, float vz) {
  TVector3 v3(vx, vy, vz);
  _vtx = v3;
}


void TCElectron::SetCharge(int c){
  _charge = c;
}

void TCElectron::SetNumberOfValidPixelHits(int n) {
  _numberOfValidPixelHits = n;
}

void TCElectron::SetNumberOfValidTrackerHits(int n) {
  _numberOfValidTrackerHits = n;
}

void TCElectron::SetNumberOfLostPixelHits(int n) {
  _numberOfLostPixelHits = n;
}

void TCElectron::SetNumberOfLostTrackerHits(int n) {
  _numberOfLostTrackerHits = n;
}

void TCElectron::SetNormalizedChi2(float n) {
  _normalizedChi2 = n;
}


void TCElectron::SetEmIso03(float e){
  _emIso03 = e;
}
void TCElectron::SetHadIso03(float h){
  _hadIso03 = h;
}
void TCElectron::SetTrkIso03(float t){
  _trkIso03 = t;
}
void TCElectron::SetEmIso04(float e){
  _emIso04 = e;
}
void TCElectron::SetHadIso04(float h){
  _hadIso04 = h;
}
void TCElectron::SetTrkIso04(float t){
  _trkIso04 = t;
}

void TCElectron::SetHadOverEm(float he){
  _hadOverEm = he;
}
void TCElectron::SetDphiSuperCluster(float dp){
  _dPhiSuperCluster = dp;
}
void TCElectron::SetDetaSuperCluster(float de){
  _dEtaSuperCluster = de;
}
void TCElectron::SetSigmaIetaIeta(float sieie){
  _sigmaIetaIeta = sieie;
}

void TCElectron::SetConversionDist(float d) {
  _convDist = d;
}

void TCElectron::SetConversionDcot(float d) {
  _convDcot = d;
}

void TCElectron::SetConversionRad(float r) {
  _convRad = r;
}

void TCElectron::SetConversionFlag(int f){
  _convFlag = f;
}

void TCElectron::SetPfSumPt(float coneSize, float f) {
  if(fabs(coneSize - 0.3) < 0.01) _pfIso_Pt03 = f;
  if(fabs(coneSize - 0.4) < 0.01) _pfIso_Pt04 = f;
  if(fabs(coneSize - 0.5) < 0.01) _pfIso_Pt05 = f;
}

void TCElectron::SetPfEGamma(float coneSize, float f) {
  if(fabs(coneSize - 0.3) < 0.01) _pfIso_Gamma03 = f;
  if(fabs(coneSize - 0.4) < 0.01) _pfIso_Gamma04 = f;
  if(fabs(coneSize - 0.5) < 0.01) _pfIso_Gamma05 = f;
}

void TCElectron::SetPfENeutral(float coneSize, float f) {
  if(fabs(coneSize - 0.3) < 0.01) _pfIso_Neutral03 = f;
  if(fabs(coneSize - 0.4) < 0.01) _pfIso_Neutral04 = f;
  if(fabs(coneSize - 0.5) < 0.01) _pfIso_Neutral05 = f;
}

void TCElectron::SetIsEB(bool b) {
  _isEB = b;
}

void TCElectron::SetIsEE(bool b) {
  _isEE = b;
}

void TCElectron::SetIsInGap(bool b) {
  _isInGap = b;
}

void TCElectron::SetCutLevel(int cut, int lvl){
  if(lvl==95){
    _cut95 = cut;
  }else if(lvl==90) {
    _cut90 = cut;
  }else if(lvl==85) {
    _cut85 = cut;
  }else if(lvl==80) {
    _cut80 = cut;
  }else if(lvl==70) {
    _cut70 = cut;
  }else if(lvl==60) {
    _cut60 = cut;
  }
}

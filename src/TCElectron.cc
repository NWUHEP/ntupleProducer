#include "TCElectron.h"
#include<iostream>

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

float TCElectron::pt() const {
   return _p4.Pt();
}

TVector3 TCElectron::Vtx() const {
   return _vtx;
}

float TCElectron::ptError() const {
  return _ptError;
}

float TCElectron::Et() const {
  return _p4.Et();
}


float TCElectron::eta() const {
 return _p4.Eta();
}

float TCElectron::phi() const {
  return _p4.Phi();
}

int TCElectron::charge() const {
   return _charge;
}

int TCElectron::numberOfValidPixelHits() const {
  return _numberOfValidPixelHits;
}

int TCElectron::numberOfValidTrackerHits() const {
  return _numberOfValidTrackerHits;
}

int TCElectron::numberOfLostPixelHits() const {
  return _numberOfLostPixelHits;
}

int TCElectron::numberOfLostTrackerHits() const {
  return _numberOfLostTrackerHits;
}

float TCElectron::normalizedChi2() const {
  return _normalizedChi2;
}


float TCElectron::emIso() const {
   return _emIso;
}

float TCElectron::hadIso() const {
   return _hadIso;
}

float TCElectron::trkIso() const {
   return _trkIso;
}

float TCElectron::pfRelIso(float coneSize) const {
  float relIso = 0;
  if (fabs(coneSize - 0.3) < 0.01)
    relIso = (_pfIso_Pt03 + _pfIso_Gamma03 + _pfIso_Neutral03) / _p4.Pt();
  if (fabs(coneSize - 0.4) < 0.01)
    relIso = (_pfIso_Pt04 + _pfIso_Gamma04 + _pfIso_Neutral04) / _p4.Pt();
  if (fabs(coneSize - 0.5) < 0.01)
    relIso = (_pfIso_Pt05 + _pfIso_Gamma05 + _pfIso_Neutral05) / _p4.Pt();
  return relIso;
}

float TCElectron::pfSumPt(float coneSize) const {
  float sumPt = 0;
  if (fabs(coneSize - 0.3) < 0.01) sumPt = _pfIso_Pt03;
  if (fabs(coneSize - 0.4) < 0.01) sumPt = _pfIso_Pt04;
  if (fabs(coneSize - 0.5) < 0.01) sumPt = _pfIso_Pt05;
  return sumPt;
}

float TCElectron::pfENeutral(float coneSize) const {
  float neutral = 0;
  if (fabs(coneSize - 0.3) < 0.01) neutral = _pfIso_Neutral03;
  if (fabs(coneSize - 0.4) < 0.01) neutral = _pfIso_Neutral04;
  if (fabs(coneSize - 0.5) < 0.01) neutral = _pfIso_Neutral05;
  return neutral;
}

float TCElectron::pfEGamma(float coneSize) const {
  float gamma = 0;
  if (fabs(coneSize - 0.3) < 0.01) gamma = _pfIso_Gamma03;
  if (fabs(coneSize - 0.4) < 0.01) gamma = _pfIso_Gamma04;
  if (fabs(coneSize - 0.5) < 0.01) gamma = _pfIso_Gamma05;
  return gamma;
}

bool TCElectron::isEB() const {
  return _isEB;
}

bool TCElectron::isEE() const {
  return _isEE;
}

bool TCElectron::isInGap() const {
  return _isInGap;
}


float TCElectron::hadOverEm() const {
  return _hadOverEm;
}
float TCElectron::dPhiSuperCluster() const {
  return _dPhiSuperCluster;
}
float TCElectron::dEtaSuperCluster() const {
  return _dEtaSuperCluster;
}
float TCElectron::sigmaIetaIeta() const {
  return _sigmaIetaIeta;
}

int TCElectron::conversionFlag() const {
  return _convFlag;
}

float TCElectron::conversionDist() const {
  return _convDist;
}

float TCElectron::conversionDcot() const {
  return _convDcot;
}

float TCElectron::conversionRad() const {
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

void TCElectron::SetnumberOfValidPixelHits(int n) {
  _numberOfValidPixelHits = n;
}

void TCElectron::SetnumberOfValidTrackerHits(int n) {
  _numberOfValidTrackerHits = n;
}

void TCElectron::SetnumberOfLostPixelHits(int n) {
  _numberOfLostPixelHits = n;
}

void TCElectron::SetnumberOfLostTrackerHits(int n) {
  _numberOfLostTrackerHits = n;
}

void TCElectron::SetnormalizedChi2(float n) {
  _normalizedChi2 = n;
}


void TCElectron::SetEMIso(float e){
  _emIso = e;
}
void TCElectron::SetHADIso(float h){
  _hadIso = h;
}
void TCElectron::SetTRKIso(float t){
  _trkIso = t;
}
void TCElectron::SetHadOverEm(float he){
  _hadOverEm = he;
}
void TCElectron::SetDPhiSuperCluster(float dp){
  _dPhiSuperCluster = dp;
}
void TCElectron::SetDEtaSuperCluster(float de){
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

void TCElectron::SetPFSumPt(float coneSize, float f) {
  if(fabs(coneSize - 0.3) < 0.01) _pfIso_Pt03 = f;
  if(fabs(coneSize - 0.4) < 0.01) _pfIso_Pt04 = f;
  if(fabs(coneSize - 0.5) < 0.01) _pfIso_Pt05 = f;
}

void TCElectron::SetPFEGamma(float coneSize, float f) {
  if(fabs(coneSize - 0.3) < 0.01) _pfIso_Gamma03 = f;
  if(fabs(coneSize - 0.4) < 0.01) _pfIso_Gamma04 = f;
  if(fabs(coneSize - 0.5) < 0.01) _pfIso_Gamma05 = f;
}

void TCElectron::SetPFENeutral(float coneSize, float f) {
  if(fabs(coneSize - 0.3) < 0.01) _pfIso_Neutral03 = f;
  if(fabs(coneSize - 0.4) < 0.01) _pfIso_Neutral04 = f;
  if(fabs(coneSize - 0.5) < 0.01) _pfIso_Neutral05 = f;
}

void TCElectron::SetisEB(bool b) {
  _isEB = b;
}

void TCElectron::SetisEE(bool b) {
  _isEE = b;
}

void TCElectron::SetisInGap(bool b) {
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

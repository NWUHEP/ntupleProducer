#include "../interface/TCEGamma.h"
#include "TCEGammaLinkDef.h"
#include <iostream>

TCEGamma::TCEGamma() {}
TCEGamma::~TCEGamma() {}

// "get" methods -------------------------------------

std::vector<TCEGamma::CrystalInfo> TCEGamma::GetCrystalVect() const { return _crysVect; }
int   TCEGamma::GetNCrystals() const { return _nCrystals;}

//float TCEGamma::MvaID() const { return _mvaID; } 
//float TCEGamma::EnergyRegression() const { return _regEne; } 
//float TCEGamma::EnergyRegressionErr() const { return _regErr; } 

float TCEGamma::R9() const { return _r9; } 


//float TCEGamma::NormalizedChi2() const {
//  return _normalizedChi2;
//}

//bool TCEGamma::ConversionVeto() const {
//    return _convVeto;
//}

//short TCEGamma::ConversionMissHits() const {
//    return _convMissHits;
//}

float TCEGamma::SCEta() const {
    return _scEta;
}

bool TCEGamma::IsEB() const {
  return _isEB;
}

bool TCEGamma::IsEE() const {
  return _isEE;
}

bool TCEGamma::IsInGap() const {
  return _isInGap;
}

float TCEGamma::HadOverEm() const {
  return _hadOverEm;
}
float TCEGamma::DphiSuperCluster() const {
  return _dPhiSuperCluster;
}
float TCEGamma::DetaSuperCluster() const {
  return _dEtaSuperCluster;
}
float TCEGamma::SigmaIEtaIEta() const {
  return _sigmaIetaIeta;
}

float TCEGamma::EOverP() const {
    return _eOverP;
}

float TCEGamma::FBrem() const {
    return _fBrem;
}

/*
bool TCEGamma::PassID(int lvl) const { 
  unsigned c = CutLevel(lvl);
  if (c & 0x01) return true;
  else return false;
}   

bool TCEGamma::PassIsolation(int lvl) const {
  unsigned c = CutLevel(lvl);
  if (c & 0x02) return true;
  else return false;
}

bool TCEGamma::PassConversion(int lvl) const {
  unsigned c = CutLevel(lvl);
  if (c & 0x04) return true;
  else return false;
}
*/

//------------------------------------------------
// "set" methods ---------------------------------------------
//------------------------------------------------------------------------



void TCEGamma::AddCrystal(TCEGamma::CrystalInfo crys) {_crysVect.push_back(crys);}
void TCEGamma::SetNCrystals(int n){ _nCrystals = n;}



void TCEGamma::SetR9(float r){ _r9 = r; } 

//void TCEGamma::SetEnergyRegression(float e){ _regEne = e; } 
//void TCEGamma::SetEnergyRegressionErr(float e){ _regErr = e; } 

void TCEGamma::SetHadOverEm(float he){
  _hadOverEm = he;
}
void TCEGamma::SetDphiSuperCluster(float dp){
  _dPhiSuperCluster = dp;
}
void TCEGamma::SetDetaSuperCluster(float de){
  _dEtaSuperCluster = de;
}

void TCEGamma::SetSigmaIEtaIEta(float sieie){
  _sigmaIetaIeta = sieie;
}

void TCEGamma::SetEOverP(float e)
{
    _eOverP = e;
}

void TCEGamma::SetFBrem(float fb)
{
    _fBrem = fb;
}

void TCEGamma::SetSCEta(float e)
{
    _scEta = e;
}

//void TCEGamma::SetConversionVeto(bool v) {
//  _convVeto = v;
//}

//void TCEGamma::SetConversionMissHits(short m) {
//  _convMissHits = m;
//}

void TCEGamma::SetIsEB(bool b) {
  _isEB = b;
}

void TCEGamma::SetIsEE(bool b) {
  _isEE = b;
}

void TCEGamma::SetIsInGap(bool b) {
  _isInGap = b;
}

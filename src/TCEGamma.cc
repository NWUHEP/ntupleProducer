#include "../interface/TCEGamma.h"
#include "TCEGammaLinkDef.h"
#include <iostream>

TCEGamma::TCEGamma() {

  //Setting up default values. In case it's not going to be set in the analyzer, 
  // one can recognize that.
  _nCrystals = 0;
  
  _r9     = -99;
  _hadOverEm = -99;

  _scEta = -99;
  _scPhi = -99;
  _scDeltaPhi = -99;
  _scDeltaEta = -99;
  _scSigmaIetaIeta  = -99;
  _scSigmaIphiIphi  = -99;

  _scEtaWidth = -99;
  _scPhiWidth = -99;
  _scEnergy   = -99;
  _e1x5 = -99;
  _e2x5 = -99;
  _e5x5 = -99;
  _preShowerOverRaw = -99;

  _pfIsoCharged = -99;
  _pfIsoNeutral = -99;
  _pfIsoPhoton  = -99;

}
TCEGamma::~TCEGamma() {}

// "get" methods -------------------------------------

std::vector<TCEGamma::CrystalInfo> TCEGamma::GetCrystalVect() const { return _crysVect; }
int   TCEGamma::GetNCrystals() const { return _nCrystals;}

//float TCEGamma::MvaID() const { return _mvaID; } 
//float TCEGamma::EnergyRegression() const { return _regEne; } 
//float TCEGamma::EnergyRegressionErr() const { return _regErr; } 

float TCEGamma::R9() const { return _r9; } 


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



//Super cluster getters 
float TCEGamma::SCEta() const {
  return _scEta;
}

float TCEGamma::SCPhi() const {
  return _scPhi;
}

float TCEGamma::SCDeltaEta() const {
  return _scDeltaEta;
}
float TCEGamma::SCDeltaPhi() const {
  return _scDeltaPhi;
}

float TCEGamma::SigmaIEtaIEta() const {
  return _scSigmaIetaIeta;
}
float TCEGamma::SigmaIPhiIPhi() const {
  return _scSigmaIphiIphi;
}

float TCEGamma::SCEtaWidth() const {
  return _scEtaWidth;
}
float TCEGamma::SCPhiWidth() const {
  return _scPhiWidth;
}

float TCEGamma::PreShowerOverRaw() const {
  return _preShowerOverRaw;
}


float TCEGamma::E1x5() const {
  return _e1x5;
}
float TCEGamma::E2x5() const {
  return _e2x5;
}
float TCEGamma::E5x5() const {
  return _e5x5;
}


float TCEGamma::SCEnergy() const {
  return _scEnergy;
}

float TCEGamma::PfIsoCharged() const {
  return _pfIsoCharged;
}
float TCEGamma::PfIsoNeutral() const {
  return _pfIsoNeutral;
}
float TCEGamma::PfIsoPhoton() const {
  return _pfIsoPhoton;
}

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

//Supercluster setters
void TCEGamma::SetSCEta(float e){
  _scEta = e;
}
void TCEGamma::SetSCPhi(float p){
  _scPhi = p;
}

void TCEGamma::SetSCDeltaEta(float de){
  _scDeltaEta = de;
}
void TCEGamma::SetSCDeltaPhi(float dp){
  _scDeltaPhi = dp;
}

void TCEGamma::SetSigmaIEtaIEta(float sieie){
  _scSigmaIetaIeta = sieie;
}

void TCEGamma::SetSigmaIPhiIPhi(float sipip){
  _scSigmaIphiIphi = sipip;
}

void TCEGamma::SetSCEtaWidth(float w){
  _scEtaWidth = w;
}
void TCEGamma::SetSCPhiWidth(float w){
  _scPhiWidth = w;
}

void TCEGamma::SetSCEnergy(float e){
  _scEnergy = e;
}

void TCEGamma::SetPreShowerOverRaw(float p){
  _preShowerOverRaw = p;
}

void TCEGamma::SetE1x5(float e){
  _e1x5 = e;
}
void TCEGamma::SetE2x5(float e){
  _e2x5 = e;
}
void TCEGamma::SetE5x5(float e){
  _e5x5 = e;
}


void TCEGamma::SetIsEB(bool b) {
  _isEB = b;
}

void TCEGamma::SetIsEE(bool b) {
  _isEE = b;
}

void TCEGamma::SetIsInGap(bool b) {
  _isInGap = b;
}


void TCEGamma::SetPfIsoCharged(float f) {
  _pfIsoCharged = f;
}
void TCEGamma::SetPfIsoNeutral(float f) {
  _pfIsoNeutral = f;
}
void TCEGamma::SetPfIsoPhoton(float f) {
  _pfIsoPhoton = f;
}

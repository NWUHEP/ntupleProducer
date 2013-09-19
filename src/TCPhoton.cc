#include "../interface/TCPhoton.h"
#include "TCPhotonLinkDef.h"
#include <iostream>


//TCPhoton::TCPhoton() { }


TCPhoton::TCPhoton() {}

/*
TCPhoton::TCPhoton(const TCPhoton& other){
  _crysArray = new CrystalInfo[100];
  for(int i=0; i < 100; i++){
    TCPhoton:: SetCrystal(i, other.GetCrystalArray()[i]);
  }
}

TCPhoton& TCPhoton::operator=(const TCPhoton& other){
  TCPhoton temp(other);
  swap(_crysArray,temp._crysArray);
  return *this;
}
*/

//TCPhoton::~TCPhoton() { delete[] _crysArray; }
TCPhoton::~TCPhoton() {}

// "get" methods -------------------------------------

std::vector<TCPhoton::CrystalInfo> TCPhoton::GetCrystalVect() const { return _crysVect; }
int   TCPhoton::GetNCrystals() const { return _nCrystals;}
float TCPhoton::NormChi2() const { return _normChi2; }
float TCPhoton::HadOverEm() const { return _hadOverEm; } 
float TCPhoton::SigmaIEtaIEta() const { return _sigmaIEtaIEta; } 
float TCPhoton::R9() const { return _r9; } 
float TCPhoton::SigmaIPhiIPhi() const { return _sigmaIPhiIPhi; } 
float TCPhoton::E2OverE9() const { return _e2OverE9; } 
bool  TCPhoton::TrackVeto() const { return _trackVeto; }

float TCPhoton::SCDPhi() const { return _SCdPhi; }
float TCPhoton::SCDEta() const { return _SCdEta; } 
float TCPhoton::SCEnergy() const { return _SCenergy; }
float TCPhoton::SCEta() const { return _SCeta; }
float TCPhoton::SCPhi() const { return _SCphi; }

bool  TCPhoton::ConversionVeto() const { return _convVeto; }

// "set" methods ---------------------------------------------


void TCPhoton::AddCrystal(TCPhoton::CrystalInfo crys) {_crysVect.push_back(crys);}
void TCPhoton::SetNCrystals(int n){ _nCrystals = n;}


void TCPhoton::SetNormChi2(float c){ _normChi2 = c; } 
void TCPhoton::SetHadOverEm(float h){ _hadOverEm = h; } 
void TCPhoton::SetSigmaIEtaIEta(float s){ _sigmaIEtaIEta = s; } 
void TCPhoton::SetSigmaIPhiIPhi(float s) { _sigmaIPhiIPhi = s; } 
void TCPhoton::SetR9(float r){ _r9 = r; } 
void TCPhoton::SetE2OverE9(float e) { _e2OverE9 = e; } 
void TCPhoton::SetTrackVeto(bool t) { _trackVeto = t; } 

void TCPhoton::SetSCDPhi(float d){ _SCdPhi = d; } 
void TCPhoton::SetSCDEta(float d){ _SCdEta = d; } 
void TCPhoton::SetSCEta(float n) { _SCeta = n; }
void TCPhoton::SetSCPhi(float p) { _SCphi = p; }
void TCPhoton::SetSCEnergy(float e) { _SCenergy = e; }

void TCPhoton::SetConversionVeto(bool v) { _convVeto = v; }


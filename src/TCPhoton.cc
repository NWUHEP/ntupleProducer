#include "../interface/TCPhoton.h"
#include "../interface/TCPhotonLinkDef.h"
#include <iostream>


//TCPhoton::TCPhoton() { }


TCPhoton::TCPhoton() {
  _crysArray = new CrystalInfo[100];

  CrystalInfo crystal;

  crystal.rawId = -99;
  crystal.ieta= -99;
  crystal.iphi= -99;
  crystal.ix= -99;
  crystal.iy= -99;
  crystal.energy= -99;
  crystal.time= -99;
  crystal.timeErr= -99;
  crystal.recoFlag= -99;

  for(int i=0; i < 100; i++){
    TCPhoton:: SetCrystal(i, crystal);
  }
  

}


TCPhoton::~TCPhoton() { delete[] _crysArray; }

// "get" methods -------------------------------------

TCPhoton::CrystalInfo* TCPhoton::GetCrystalArray() const { return _crysArray; }
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


void TCPhoton::SetCrystal(int i , CrystalInfo crys) {
  _crysArray[i].rawId = crys.rawId;
  _crysArray[i].ieta =crys.ieta;
  _crysArray[i].iphi =crys.iphi;
    _crysArray[i].ix =crys.ix;
  _crysArray[i].iy =crys.iy;
  _crysArray[i].energy =crys.energy;
  _crysArray[i].time =crys.time;
  _crysArray[i].timeErr =crys.timeErr;
  _crysArray[i].recoFlag =crys.recoFlag;
}
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


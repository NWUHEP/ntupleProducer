#include "../interface/TCElectron.h"
#include "TCElectronLinkDef.h"
#include <iostream>

TCElectron::TCElectron() {
  _normChi2gsf = -99;
  _normChi2kf  = -99;
  
  _ptError     = -99;
  

}

TCElectron::~TCElectron() {
}



// "get" methods -------------------------------------
float TCElectron::NormalizedChi2() const { 
  return _normChi2gsf; 
}
float TCElectron::NormalizedChi2Gsf() const { 
  return _normChi2gsf; 
}
float TCElectron::NormalizedChi2Kf()  const { 
  return _normChi2kf; 
}

float TCElectron::InverseEnergyMomentumDiff() const{
  return _inverseEnergyMomentumDiff;
}


float TCElectron::MvaID() const { 
  return _mvaID; 
} 
float TCElectron::EnergyRegression() const { 
  return _regEne; 
}
float TCElectron::EnergyRegressionErr() const { 
  return _regErr; 
} 

float TCElectron::PtError() const {
  return _ptError;
}

int TCElectron::TrackerLayersWithMeasurement()  const {
  return _trackerLayersWithMeasurement; 
}

int TCElectron::NumberOfValidHits() const {
  return _numberOfValidHits;
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


bool TCElectron::ConversionVeto() const {
    return _convVeto;
}

short TCElectron::ConversionMissHits() const {
    return _convMissHits;
}


float TCElectron::FBrem() const {
    return _fBrem;
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
}   

bool TCElectron::PassIsolation(int lvl) const {
  unsigned c = CutLevel(lvl);
  if (c & 0x02) return true;
  else return false;
}

bool TCElectron::PassConversion(int lvl) const {
  unsigned c = CutLevel(lvl);
  if (c & 0x04) return true;
  else return false;
}

TLorentzVector TCElectron::RegressionMomCombP4() const {
  return _regressionMomCombP4;
}

float TCElectron::EffArea() const {
  return _effArea;
}
//------------------------------------------------
// "set" methods ---------------------------------------------
//------------------------------------------------------------------------

void TCElectron::SetNormalizedChi2Gsf(float c){ 
  _normChi2gsf = c; 
} 
void TCElectron::SetNormalizedChi2Kf(float c) { 
  _normChi2kf  = c; 
} 

void TCElectron::SetInverseEnergyMomentumDiff(float d){
  _inverseEnergyMomentumDiff = d;
}

void TCElectron::SetMvaID(float m){ 
  _mvaID = m; 
}
void TCElectron::SetEnergyRegression(float e){ 
  _regEne = e; 
}
void TCElectron::SetEnergyRegressionErr(float e){ 
  _regErr = e; 
} 

void TCElectron::SetPtError(float e) {
  _ptError = e;
}

void TCElectron::SetTrackerLayersWithMeasurement(int t) { 
  _trackerLayersWithMeasurement = t; 
}

void TCElectron::SetNumberOfValidHits(int n) {
  _numberOfValidHits = n;
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

void TCElectron::SetFBrem(float fb)
{
    _fBrem = fb;
}


void TCElectron::SetConversionVeto(bool v) {
  _convVeto = v;
}

void TCElectron::SetConversionMissHits(short m) {
  _convMissHits = m;
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

void TCElectron::SetRegressionMomCombP4(TLorentzVector tmpP4){
  _regressionMomCombP4 = tmpP4;
}

void TCElectron::SetEffArea(float a){
  _effArea = a;
}


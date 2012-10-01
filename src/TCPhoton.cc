#include "../interface/TCPhoton.h"
#include <iostream>

TCPhoton::TCPhoton() { }

TCPhoton::~TCPhoton() { }

// "get" methods -------------------------------------

float TCPhoton::NormChi2() const { return _normChi2; }
float TCPhoton::HadOverEm() const { return _hadOverEm; } 
float TCPhoton::SigmaIEtaIEta() const { return _sigmaIEtaIEta; } 
float TCPhoton::R9() const { return _r9; } 
float TCPhoton::SigmaIPhiIPhi() const { return _sigmaIPhiIPhi; } 
float TCPhoton::E2OverE9() const { return _e2OverE9; } 
bool  TCPhoton::TrackVeto() const { return _trackVeto; }

float TCPhoton::DPhiSC() const { return _dPhiSC; }
float TCPhoton::DEtaSC() const { return _dEtaSC; } 
float TCPhoton::EnergySC() const { return _energySC; }
float TCPhoton::EtaSC() const { return _etaSC; }

bool  TCPhoton::ConversionVeto() const { return _convVeto; }

// "set" methods ---------------------------------------------

void TCPhoton::SetNormChi2(float c){ _normChi2 = c; } 
void TCPhoton::SetHadOverEm(float h){ _hadOverEm = h; } 
void TCPhoton::SetSigmaIEtaIEta(float s){ _sigmaIEtaIEta = s; } 
void TCPhoton::SetSigmaIPhiIPhi(float s) { _sigmaIPhiIPhi = s; } 
void TCPhoton::SetR9(float r){ _r9 = r; } 
void TCPhoton::SetE2OverE9(float e) { _e2OverE9 = e; } 
void TCPhoton::SetTrackVeto(bool t) { _trackVeto = t; } 

void TCPhoton::SetDPhiSC(float d){ _dPhiSC = d; } 
void TCPhoton::SetDEtaSC(float d){ _dEtaSC = d; } 
void TCPhoton::SetEtaSC(float n) { _etaSC = n; }
void TCPhoton::SetEnergySC(float e) { _energySC = e; }

void TCPhoton::SetConversionVeto(bool v) { _convVeto = v; }

#include "../interface/TCPhoton.h"
#include "TCPhotonLinkDef.h"
#include <iostream>

TCPhoton::TCPhoton() {

}

TCPhoton::~TCPhoton() {}

// "get" methods -------------------------------------
//float TCPhoton::E2OverE9() const { return _e2OverE9; } 
bool  TCPhoton::TrackVeto() const { 
  return _trackVeto; 
}

bool  TCPhoton::ConversionVeto() const { 
  return _convVeto; 
}

int TCPhoton::NTrkSolidConeDR03() const {
  return _nTrkSolidConeDR03;
}

// "set" methods ---------------------------------------------

//void TCPhoton::SetE2OverE9(float e) { _e2OverE9 = e; } 
void TCPhoton::SetTrackVeto(bool t) { 
  _trackVeto = t; 
} 
void TCPhoton::SetConversionVeto(bool v) { 
  _convVeto = v; 
}
void TCPhoton::SetNTrkSolidConeDR03(int n){
  _nTrkSolidConeDR03 = n;
}

#include "../interface/TCTrack.h"
#include "TCTrackLinkDef.h"

TCTrack::TCTrack() {
  _normChi2 = -99;
}

TCTrack::~TCTrack() {}

float TCTrack::NormalizedChi2() const {
  return _normChi2;
}

void TCTrack::SetNormalizedChi2(float c){
  _normChi2 = c;
} 

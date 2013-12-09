#include "../interface/TCTrack.h"
#include "TCTrackLinkDef.h"

TCTrack::TCTrack() {
  _normChi2 = -99;
}

TCTrack::~TCTrack() {}

float TCTrack::NormalizedChi2() const {
  return _normChi2;
}
float TCTrack::PtError() const {
  return _ptError;
}

void TCTrack::SetNormalizedChi2(float c){
  _normChi2 = c;
}

void TCTrack::SetPtError(float e){
  _ptError = e;
}

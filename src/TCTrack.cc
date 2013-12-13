#include "../interface/TCTrack.h"
#include "TCTrackLinkDef.h"

TCTrack::TCTrack():
  _normChi2(-99),_ptError(-99)
{}

TCTrack::~TCTrack() {}

TCTrack::ConversionInfo TCTrack::GetConversionInfo() const {
  return _convInfo;
}

float TCTrack::NormalizedChi2() const {
  return _normChi2;
}
float TCTrack::PtError() const {
  return _ptError;
}

void TCTrack::SetConversionInfo(TCTrack::ConversionInfo i){
  _convInfo = i;
}

void TCTrack::SetNormalizedChi2(float c){
  _normChi2 = c;
}

void TCTrack::SetPtError(float e){
  _ptError = e;
}

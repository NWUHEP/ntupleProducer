#ifndef _TCTRACK_H
#define	_TCTRACK_H

#include "TCPhysObject.h"

class TCTrack : public TCPhysObject {

 public:
  struct ConversionInfo{
    bool isValid;
    int nHitsMax;
    double vtxProb;
    double lxyPV;
    double lxyBS;
    
    ConversionInfo() : isValid(false), nHitsMax(-1), vtxProb(-1), lxyPV(-99), lxyBS(-99) {}
  };


 private:
  float _normChi2;
  float _ptError;
  TCTrack::ConversionInfo _convInfo;

 public:
  TCTrack();
  virtual ~TCTrack();

  TCTrack::ConversionInfo GetConversionInfo() const;
  float NormalizedChi2 () const;
  float PtError () const;

  void SetConversionInfo(TCTrack::ConversionInfo);
  void SetNormalizedChi2 (float);
  void SetPtError (float);

  ClassDef(TCTrack, 1);

};

#endif	/* _TCTRACK_H */



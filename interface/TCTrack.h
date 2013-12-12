#ifndef _TCTRACK_H
#define	_TCTRACK_H

#include "TCPhysObject.h"

class TCTrack : public TCPhysObject {
public 
 private:
  float _normChi2;
  float _ptError;

 public:
  TCTrack();
  virtual ~TCTrack();

  float NormalizedChi2 () const;
  float PtError () const;

  void SetNormalizedChi2 (float);
  void SetPtError (float);

  ClassDef(TCTrack, 1);

};

#endif	/* _TCTRACK_H */



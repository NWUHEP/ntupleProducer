#ifndef _TCTRACK_H
#define	_TCTRACK_H

#include "TCPhysObject.h"

class TCTrack : public TCPhysObject {
 private:
  float _normChi2;
  
 public:
  TCTrack();
  virtual ~TCTrack();
  
  float NormalizedChi2 () const;


  void SetNormalizedChi2 (float);

  ClassDef(TCTrack, 1);

};

#endif	/* _TCTRACK_H */



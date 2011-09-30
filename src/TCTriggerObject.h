#ifndef _TCTriggerObject_H
#define _TCTriggerObject_H

#include "TObject.h"
#include "TLorentzVector.h"

class TCTriggerObject : public TObject {
 private:
  TLorentzVector _p4;
  int _id;
  
 public:
  TCTriggerObject();
  virtual ~TCTriggerObject();
  
  void SetId(int i);
  void SetP4(double px, double py, double pz, double energy);
  
  TLorentzVector P4();
  int Id();
  
  ClassDef(TCTriggerObject, 1);
};

#endif  /* _TCTriggerObject_H  */


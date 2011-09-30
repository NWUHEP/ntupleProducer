#include "TCTriggerObject.h"

TCTriggerObject::TCTriggerObject() {
}

TCTriggerObject::~TCTriggerObject() {
}

void TCTriggerObject::SetId(int i) {
  _id = i;
}

void TCTriggerObject::SetP4(double px, double py, double pz, double energy) {
  TLorentzVector blah(px, py, pz, energy);
  _p4 = blah;
}

TLorentzVector TCTriggerObject::P4() {
  return _p4;
}

int TCTriggerObject::Id() {
  return _id;
}

#include "TCTriggerObject.h"

TCTriggerObject::TCTriggerObject() {
}

TCTriggerObject::~TCTriggerObject() {
}

void TCTriggerObject::setId(int i) {
       id = i;
}

void TCTriggerObject::setP4(double px, double py, double pz, double energy) {
       TLorentzVector blah(px, py, pz, energy);
          p4 = blah;
}

TLorentzVector TCTriggerObject::P4() {
       return p4;
}

int TCTriggerObject::getId() {
       return id;
}

#include "../interface/TCGenParticle.h"
#include "TCGenParticleLinkDef.h"

TCGenParticle::TCGenParticle() {
  _mother = 0;
  _PDGID = 0;
  _status = 0;
  _isParton = 0;
}

TCGenParticle::~TCGenParticle() {
}

// "get" methods -------------------------------------
/*
Needs more validation. May not work in the producer code, too many references
TCGenParticle* TCGenParticle::PrimaryAncestor(){
  TCGenParticle *a = this;
  while (a->Mother())
    a = a->Mother();
  return a;
}
*/

TCGenParticle* TCGenParticle::Mother() {
    return _mother;
}

int TCGenParticle::GetPDGId() {
    return _PDGID;
}

unsigned TCGenParticle::GetStatus() {
    return _status;
}

bool TCGenParticle::IsParton() {
  return _isParton;
}

void TCGenParticle::SetMother(TCGenParticle* m) {
    _mother = m;
}

void TCGenParticle::SetPDGId(int p) {
    _PDGID = p;
}

void TCGenParticle::SetStatus(unsigned s)
{
    _status = s;
}
void TCGenParticle::SetIsParton(bool a) {
    _isParton = a;
}

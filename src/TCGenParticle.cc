#include "../interface/TCGenParticle.h"
#include "TCGenParticleLinkDef.h"

TCGenParticle::TCGenParticle():
  _momID(-99),
  _mother(0),
  _PDGID(-99),
  _status(-99),
  _isParton(0)
{
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

int TCGenParticle::GetPDGId() const{
    return _PDGID;
}

int TCGenParticle::MotherId() const{
    return _momID;
}

unsigned TCGenParticle::GetStatus() const{
    return _status;
}

bool TCGenParticle::IsParton() const{
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

void TCGenParticle::SetMotherId(int id) {
    _momID = id;
}

ostream& TCGenParticle::TCprint(ostream& os) const {
 return TCPhysObject::TCprint(os) << 
   " PDGId: " << GetPDGId() << " Status: " << GetStatus() << " IsParton: " << IsParton() << " MotherId: " << MotherId();
}

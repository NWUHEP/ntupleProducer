#include "../interface/TCGenParticle.h"
#include "TCGenParticleLinkDef.h"
#include <iostream>

TCGenParticle::TCGenParticle() {
  mother = 0;
}

TCGenParticle::~TCGenParticle() {
}

// "get" methods -------------------------------------

TCGenParticle* TCGenParticle::Mother() {
    return mother;
}

int TCGenParticle::GetPDGId() {
    return PDGID;
}

unsigned TCGenParticle::GetStatus() {
    return status;
}

bool TCGenParticle::IsParton() {
    return isParton_;
}


void TCGenParticle::SetMother(TCGenParticle* m) {
    mother = m;
}

void TCGenParticle::SetPDGId(int pdg_id) {
    PDGID = pdg_id;
}

void TCGenParticle::SetStatus(unsigned s)
{
    status = s;
}
void TCGenParticle::SetIsParton(bool a) {
    isParton_=a;
}

#include "../interface/TCGenParticle.h"
#include <iostream>

TCGenParticle::TCGenParticle() {
}

TCGenParticle::~TCGenParticle() {
}

// "get" methods -------------------------------------

int TCGenParticle::Mother() {
    return mother;
}

int TCGenParticle::Grandmother() {
    return grandmother;
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

//std::vector<int> TCGenParticle::GetDaughters() {
//   return daughters;
//}

// "set" methods ---------------------------------------------

//void TCGenParticle::AddDaughter(int d) {
//   daughters.push_back(d);
//}

void TCGenParticle::SetMother(int m) {
    mother = m;
}

void TCGenParticle::SetGrandmother(int g) {
    grandmother = g;
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

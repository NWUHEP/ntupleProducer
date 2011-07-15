#include "TCGenParticle.h"
#include <iostream>

TCGenParticle::TCGenParticle() {
}

TCGenParticle::~TCGenParticle() {
}

// "get" methods -------------------------------------

TVector3 TCGenParticle::Position() const {
	   return _position;
}

TLorentzVector TCGenParticle::P4() const {
	   return _p4;
}

TVector2 TCGenParticle::P2() const {
	   TVector2 v2(_p4.Px(), _p4.Py());
		   return v2;
}

float TCGenParticle::Et() const {
	   return _p4.Et();
}

float TCGenParticle::Pt() const {
	   return _p4.Pt();
}

int TCGenParticle::Charge() const{
	   return charge;
}

int TCGenParticle::Mother() {
	   return mother;
}

int TCGenParticle::GetPDGId() {
	   return PDGID;
}

//std::vector<int> TCGenParticle::GetDaughters() {
//   return daughters;
//}

// "set" methods ---------------------------------------------

void TCGenParticle::SetPosition(float x, float y, float z) {
	   TVector3 p(x, y, z);
		   _position = p;
}

void TCGenParticle::SetP4(TLorentzVector p4) {
	   _p4 = p4;
}

void TCGenParticle::SetP4(float px, float py, float pz, float e) {
	   TLorentzVector p4(px, py, pz, e);
		   _p4 = p4;
}

void TCGenParticle::SetCharge(int c) {
	   charge = c;
}

//void TCGenParticle::AddDaughter(int d) {
//   daughters.push_back(d);
//}

void TCGenParticle::SetMother(int m) {
	   mother = m;
}

void TCGenParticle::SetPDGId(int pdg_id) {
	   PDGID = pdg_id;
}

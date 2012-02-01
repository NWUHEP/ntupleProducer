#include "Higgs/ntupleProducer/interface/TCTau.h"
#include <iostream>

TCTau::TCTau() {
}

TCTau::~TCTau() {
}

TVector3 TCTau::PositionFromTrack() const {
	   return _positionFromLeadTrack;
}

TVector3 TCTau::PositionFromTau() const {
	   return _positionFromTauObject;
}

TLorentzVector TCTau::P4() const {
	   return _p4;
}

TVector2 TCTau::P2() const {
	   TVector2 v2(_p4.Px(), _p4.Py());
		   return v2;
}

float TCTau::Et() const {
	   return _p4.Et();
}

float TCTau::Pt() const {
	   return _p4.Pt();
}



void TCTau::SetPositionFromTau(float x, float y, float z) {
	   TVector3 p(x, y, z);
		   _positionFromTauObject = p;
}

void TCTau::SetPositionFromTrack(float x, float y, float z) {
	   TVector3 p(x, y, z);
		   _positionFromLeadTrack = p;
}

void TCTau::SetP4(TLorentzVector p4) {
	   _p4 = p4;
}

void TCTau::SetP4(float px, float py, float pz, float e) {
	   TLorentzVector p4(px, py, pz, e);
		   _p4 = p4;
}


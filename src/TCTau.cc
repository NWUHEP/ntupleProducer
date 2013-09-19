#include "../interface/TCTau.h"
#include "TCTauLinkDef.h"
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

void TCTau::SetPositionFromTau(float x, float y, float z) {
	   TVector3 p(x, y, z);
		   _positionFromTauObject = p;
}

void TCTau::SetPositionFromTrack(float x, float y, float z) {
	   TVector3 p(x, y, z);
		   _positionFromLeadTrack = p;
}


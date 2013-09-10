/* 
 * File:   TCGenJet.cc
 * Author: Nate O.
 * 
 * Created on Nov 29, 2010, 1:08 PM
 */

#include "../interface/TCGenJet.h"
#include "../interface/TCGenJetLinkDef.h"

TCGenJet::TCGenJet() {}

TCGenJet::~TCGenJet() {}

// "get" methods -------------------------------------

float TCGenJet::HadEnergy() const {
   return _hadEnergy;
}

float TCGenJet::EmEnergy() const {
   return _emEnergy;
}

float TCGenJet::InvEnergy() const {
   return _invEnergy;
}

float TCGenJet::AuxEnergy() const {
   return _auxEnergy;
}

unsigned int TCGenJet::NumConstit() const {
   return _numConstit;
}

unsigned int TCGenJet::NumChPart() const {
   return _numChPart;
}

TLorentzVector TCGenJet::ProgenitorP4() const {
   return _progenitorP4;
}

int TCGenJet::JetFlavor() const {
   return _jetFlavor;
}

// "set" methods ---------------------------------------------

void TCGenJet::SetProgenitorP4(TLorentzVector p4) {
   _progenitorP4 = p4;
}

void TCGenJet::SetHadEnergy(float h) {
   _hadEnergy = h;
}

void TCGenJet::SetEmEnergy(float e) {
   _emEnergy = e;
}

void TCGenJet::SetInvEnergy(float i) {
   _invEnergy = i;
}

void TCGenJet::SetAuxEnergy(float a) {
   _auxEnergy = a;
}

void TCGenJet::SetNumConstit(unsigned int n) {
   _numConstit = n;
}

void TCGenJet::SetNumChPart(unsigned int n) {
   _numChPart = n;
}

void TCGenJet::SetJetFlavor(int f) {
   _jetFlavor = f;
}

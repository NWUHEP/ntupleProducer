/* 
 * File:   TCGenJet.cc
 * Author: Nate O.
 * 
 * Created on Nov 29, 2010, 1:08 PM
 */

#include "TCGenJet.h"
#include<iostream>

TCGenJet::TCGenJet() {}

TCGenJet::~TCGenJet() {}

// "get" methods -------------------------------------

TLorentzVector TCGenJet::P4() const {
   return _p4;
}

TVector2 TCGenJet::P2() const {
   TVector2 v2(_p4.Px(), _p4.Py());
   return v2;
}

float TCGenJet::Et() const {
   return _p4.Et();
}

float TCGenJet::Pt() const {
   return _p4.Pt();
}

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

TVector3 TCGenJet::Vtx() const {
   return _vtx;
}

// "set" methods ---------------------------------------------

void TCGenJet::SetP4(TLorentzVector p4) {
   _p4 = p4;
}

void TCGenJet::SetP4(float px, float py, float pz, float e) {
   TLorentzVector p4(px, py, pz, e);
   _p4 = p4;
}

void TCGenJet::SetVtx(float vx, float vy, float vz) {
   TVector3 v3(vx, vy, vz);
   _vtx = v3;
}

void TCGenJet::SetProgenitorP4(TLorentzVector p4) {
   _progenitorP4 = p4;
}

//void TCGenJet::SetAssocVtx(float vx, float vy, float vz) {
//   TVector3 v3(vx, vy, vz);
//   _assocPV = v3;
//}

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

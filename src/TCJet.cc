/* 
 * File:   TCJet.cc
 * Author: Anton A.
 * 
 * Created on April 30, 2010, 2:49 PM
 */

#include "TCJet.h"
#include<iostream>

TCJet::TCJet() {
   for (int i = 0; i < 8; ++i) {
      _jetCorr[i] = 1.0;
      _jetCorrIsSet[i] = false;
   }
   _jetCorrIsSet[0] = true;
}

TCJet::~TCJet() {
}

// "get" methods -------------------------------------

TLorentzVector TCJet::P4() const {
   return _p4;
}

TVector2 TCJet::P2() const {
   TVector2 v2(_p4.Px(), _p4.Py());
   return v2;
}

float TCJet::Et() const {
   return _p4.Et();
}

float TCJet::Pt() const {
   return _p4.Pt();
}

// accessors for corrected jets (argument is level of correction)

TLorentzVector TCJet::P4(unsigned int lvl) const {
   if (lvl > 7) {
      std::cout << "\nJet correction = " << lvl << std::endl;
      std::cout << "Correction level cannot exceed 7!\n";
      std::cout << "No correction will be applied!\n";
      return _p4;
   }
   return TotalJetCorr(lvl) * _p4;
}

TVector2 TCJet::P2(unsigned int lvl) const {
   TVector2 v2(_p4.Px(), _p4.Py());
   if (lvl > 7) {
      std::cout << "\nJet correction = " << lvl << std::endl;
      std::cout << "Correction level cannot exceed 7!\n";
      std::cout << "No correction will be applied!\n";
      return v2;
   }
   return TotalJetCorr(lvl) * v2;
}

float TCJet::Et(unsigned int lvl) const {

   if (lvl > 7) {
      std::cout << "\nJet correction = " << lvl << std::endl;
      std::cout << "Correction level cannot exceed 7!\n";
      std::cout << "No correction will be applied!\n";
      return _p4.Et();
   }
   return TotalJetCorr(lvl) * _p4.Et();
}

float TCJet::Pt(unsigned int lvl) const {
   if (lvl > 7) {
      std::cout << "\nJet correction = " << lvl << std::endl;
      std::cout << "Correction level cannot exceed 7!\n";
      std::cout << "No correction will be applied!\n";
      return _p4.Pt();
   }
   return TotalJetCorr(lvl) * _p4.Pt();
}

float TCJet::ChHadFrac() const {
   return _chHadFrac;
}

float TCJet::NeuHadFrac() const {
   return _neuHadFrac;
}

float TCJet::ChEmFrac() const {
   return _chEmFrac;
}

float TCJet::NeuEmFrac() const {
   return _neuEmFrac;
}

unsigned int TCJet::NumConstit() const {
   return _numConstit;
}

unsigned int TCJet::NumChPart() const {
   return _numChPart;
}

TVector3 TCJet::Vtx() const {
   return _vtx;
}

float TCJet::VtxSumPtFrac() const {
   return _vtxSumPtFrac;
}

float TCJet::VtxSumPt() const {
   return _vtxSumPt;
}  

float TCJet::VtxTrackFrac() const {
   return _vtxTrackFrac;
}  

int TCJet::VtxNTracks() const {
   return _vtxNTracks;
}  

unsigned int TCJet::VtxSumPtIndex() const {
   return _vtxSumPtIndex;
}

unsigned int TCJet::VtxCountIndex() const {
   return _vtxCountIndex;
}

bool TCJet::JetCorrIsSet(unsigned int lvl) const {
   return _jetCorrIsSet[lvl];
}

float TCJet::JetCorr(unsigned int lvl) const {
   return _jetCorr[lvl];
}

float TCJet::TotalJetCorr(unsigned int lvl) const {
   float corr = 1.0;
   for (unsigned int i = 1; i <= lvl; ++i) {
      if (JetCorrIsSet(lvl)) corr *= JetCorr(i);
   }
   return corr;
}

float TCJet::UncertaintyJES() const {
	return _jesUncertainty;
}

// b tagging discriminators
//Track counting tag with N = 3: trackCountingHighPurBJetTags

float TCJet::BDiscrTCHP() const {
   return _bDiscrTCHP;
}

//Track counting tag with N = 2: trackCountingHighEffBJetTags

float TCJet::BDiscrTCHE() const {
   return _bDiscrTCHE;
}

//Simple secondary vertex b tag: simpleSecondaryVertexBJetTags

float TCJet::BDiscrSSVHE() const {
   return _bDiscrSSVHE;
}

float TCJet::BDiscrSSVHP() const {
   return _bDiscrSSVHP;
}

float TCJet::BDiscrJP() const {
   return _bDiscrJP;
}

float TCJet::BDiscrJBP() const {
   return _bDiscrJBP;
}

float TCJet::BDiscrCSV() const {
   return _bDiscrCSV;
}

int TCJet::JetFlavor() const {
    return _jetFlavor;
}

// "set" methods ---------------------------------------------

void TCJet::SetP4(TLorentzVector p4) {
   _p4 = p4;
}

void TCJet::SetP4(float px, float py, float pz, float e) {
   TLorentzVector p4(px, py, pz, e);
   _p4 = p4;
}

void TCJet::SetVtx(float vx, float vy, float vz) {
   TVector3 v3(vx, vy, vz);
   _vtx = v3;
}

void TCJet::SetVtxSumPtFrac(float f){
   _vtxSumPtFrac = f;
}  

void TCJet::SetVtxSumPt(float p){
   _vtxSumPt = p;
}  

void TCJet::SetVtxTrackFrac(float f){
   _vtxTrackFrac = f;
}  

void TCJet::SetVtxNTracks(int n){
   _vtxNTracks = n;
}  

void TCJet::SetVtxSumPtIndex(unsigned int i){
   _vtxSumPtIndex = i;
} 

void TCJet::SetVtxCountIndex(unsigned int i){
   _vtxCountIndex = i;
} 

void TCJet::SetChHadFrac(float c) {
   _chHadFrac = c;
}

void TCJet::SetNeuHadFrac(float n) {
   _neuHadFrac = n;
}

void TCJet::SetChEmFrac(float c) {
   _chEmFrac = c;
}

void TCJet::SetNeuEmFrac(float n) {
   _neuEmFrac = n;
}

void TCJet::SetNumConstit(unsigned int n) {
   _numConstit = n;
}

void TCJet::SetNumChPart(unsigned int n) {
   _numChPart = n;
}

void TCJet::SetJetCorr(unsigned int lvl, float corr) {

   if (lvl <= 7) {
      _jetCorr[lvl] = corr;
      _jetCorrIsSet[lvl] = true;

   } else {
      std::cout << "\nJet correction lvl = " << lvl << " is not valid!\n";
      std::cout << "No correction will be applied!\n\n";
   }
}

void TCJet::SetUncertaintyJES(float u) {
	_jesUncertainty = u;
}
// b tagging discriminators

void TCJet::SetBDiscrTCHE(float d) {
   _bDiscrTCHE = d;
}

void TCJet::SetBDiscrTCHP(float d) {
   _bDiscrTCHP = d;
}

void TCJet::SetBDiscrSSVHE(float d) {
   _bDiscrSSVHE = d;
}

void TCJet::SetBDiscrSSVHP(float d) {
   _bDiscrSSVHP = d;
}

void TCJet::SetBDiscrJP(float d) {
   _bDiscrJP = d;
}

void TCJet::SetBDiscrJBP(float d) {
   _bDiscrJBP = d;
}

void TCJet::SetBDiscrCSV(float d) {
   _bDiscrCSV = d;
}

void TCJet::SetJetFlavor(float f) {
    _jetFlavor = f;
}

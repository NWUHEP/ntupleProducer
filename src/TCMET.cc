/* 
 * File:   TCMET.cc
 * Author: Nate O. 
 * 
 * Created on December 6 2010 8:04 PM
 */

#include "TCMET.h"
#include<iostream>

TCMET::TCMET() {}

TCMET::~TCMET() {}

// "get" methods ------------------------------------

float TCMET::SumEt() const {
  return _sumEt;
}

float TCMET::Met() const {
  return _met;
}

float TCMET::Phi() const {
  return _phi;
}
float TCMET::CorrectedSumEt() const {
  return _corSumEt;
}

float TCMET::CorrectedMet() const {
  return _corMet;
}

float TCMET::CorrectedPhi() const {
  return _corPhi;
}
float TCMET::PhotonEtFraction() const {
  return _photonEtFraction;
}

float TCMET::ElectronEtFraction() const {
  return _electronEtFraction;
}

float TCMET::MuonEtFraction() const {
  return _muonEtFraction;
}

float TCMET::NeutralHadronEtFraction() const {
  return _neutralHadronEtFraction;
}

float TCMET::ChargedHadronEtFraction() const {
  return _chargedHadronEtFraction;
}

float TCMET::HFEMEtFraction() const {
  return _hfEMEtFraction;
}

float TCMET::HFHadronEtFraction() const {
  return _hfHadronEtFraction;
}

float TCMET::Significance() const {
  if (_sumEt!=0) return _met/_sumEt;
  else return -1;
}

    // "set" methods ---------

void TCMET::SetSumEt(float n) {
  _sumEt = n;
}

void TCMET::SetMet(float n) {
  _met = n;
}

void TCMET::SetPhi(float n) {
  _phi = n;
}

void TCMET::SetCorrectedSumEt(float n) {
  _corSumEt = n;
}

void TCMET::SetCorrectedMet(float n) {
  _corMet = n;
}

void TCMET::SetCorrectedPhi(float n) {
  _corPhi = n;
}
void TCMET::SetPhotonEtFraction(float n) {
  _photonEtFraction = n;
}

void TCMET::SetElectronEtFraction(float n) {
  _electronEtFraction = n;
}

void TCMET::SetMuonEtFraction(float n) {
  _muonEtFraction = n;
}

void TCMET::SetNeutralHadronEtFraction(float n) {
  _neutralHadronEtFraction = n;
}

void TCMET::SetChargedHadronEtFraction(float n) {
  _chargedHadronEtFraction = n;
}

void TCMET::SetHFEMEtFraction(float n) {
  _hfEMEtFraction = n;
}

void TCMET::SetHFHadronEtFraction(float n) {
  _hfHadronEtFraction = n;
}




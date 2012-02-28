/* 
 * File:   TCMET.cc
 * Author: Nate O. 
 * 
 * Created on December 6 2010 8:04 PM
 */

#include "Higgs/ntupleProducer/interface/TCMET.h"
#include <iostream>

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

float TCMET::MuonFraction() const {
  return _muonFraction;
}

float TCMET::NeutralHadronFraction() const {
  return _neutralHadronFraction;
}

float TCMET::NeutralEMFraction() const {
  return _neutralEMFraction;
}

float TCMET::ChargedHadronFraction() const {
  return _chargedHadronFraction;
}

float TCMET::ChargedEMFraction() const {
  return _chargedEMFraction;
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

void TCMET::SetMuonFraction(float n) {
  _muonFraction = n;
}

void TCMET::SetNeutralHadronFraction(float n) {
  _neutralHadronFraction = n;
}

void TCMET::SetNeutralEMFraction(float n) {
  _neutralEMFraction = n;
}

void TCMET::SetChargedHadronFraction(float n) {
  _chargedHadronFraction = n;
}

void TCMET::SetChargedEMFraction(float n) {
  _chargedEMFraction = n;
}


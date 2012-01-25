/* 
 * File:   TCPrimaryVtx.cc
 * Author: Anton A.
 * 
 * Created on May 21, 2010, 11:16 AM
 */

#include "Higgs/ntupleProducer/interface/TCPrimaryVtx.h"

TCPrimaryVtx::TCPrimaryVtx() {
}

TCPrimaryVtx::~TCPrimaryVtx() {
}

TVector3 TCPrimaryVtx::Position() const {
   return _position;
}

float TCPrimaryVtx::NDof() const {
   return _nDof;
}

float TCPrimaryVtx::Chi2() const {
   return _chi2;
}

bool TCPrimaryVtx::IsFake() const {
   return _isFake;
}

int TCPrimaryVtx::Ntracks() const {
   return _nTracks;
}

float TCPrimaryVtx::SumPt2Trks() const {
  return _sumPt2Trks;
}


void TCPrimaryVtx::SetPosition(float x, float y, float z) {
   TVector3 p(x, y, z);
   _position = p;
}

void TCPrimaryVtx::SetNDof(float n) {
   _nDof = n;
}

void TCPrimaryVtx::SetChi2(float chi2) {
   _chi2 = chi2;
}

void TCPrimaryVtx::SetIsFake(bool isF) {
   _isFake = isF;
}

void TCPrimaryVtx::SetNtracks(int nTrk) {
   _nTracks = nTrk;
}

void TCPrimaryVtx::SetSumPt2Trks(float sumPt2) {
  _sumPt2Trks = sumPt2;
}

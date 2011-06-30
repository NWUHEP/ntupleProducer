#include "TCMuon.h"
#include<iostream>

TCMuon::TCMuon() {

}

TCMuon::~TCMuon() {

}



// "get" methods -------------------------------------

TLorentzVector TCMuon::P4() const {
   return _p4;
}

TVector2 TCMuon::P2() const {
  TVector2 v2(_p4.Px(), _p4.Py());
  return v2;
}

float TCMuon::Pt() const {
   return _p4.Pt();
}

float TCMuon::PtError() const {
   return _ptError;
}

float TCMuon::Et() const {
  return _p4.Et();
}

float TCMuon::Eta() const {
  return _p4.Eta();
}

float TCMuon::Phi() const {
  return _p4.Phi();
}

TVector3 TCMuon::Vtx() const {
   return _vtx;
}

int TCMuon::Charge() const {
   return _charge;
}

bool TCMuon::IsTRK() const {
   return _isTRK;
}

bool TCMuon::IsGLB() const {
   return _isGLB;
}

int TCMuon::NumberOfValidPixelHits() const {
  return _numberOfValidPixelHits;
}

int TCMuon::NumberOfValidTrackerHits() const {
  return _numberOfValidTrackerHits;
}

int TCMuon::NumberOfValidMuonHits() const {
  return _numberOfValidMuonHits;
}

int TCMuon::NumberOfLostPixelHits() const {
  return _numberOfLostPixelHits;
}

int TCMuon::NumberOfLostTrackerHits() const {
  return _numberOfLostTrackerHits;
}

float TCMuon::NormalizedChi2() const {
  return _normalizedChi2;
}

int TCMuon::NumberOfMatches() const {
   return _numberOfMatches;
}

float TCMuon::CaloComp() const {
   return _caloComp;
}

float TCMuon::SegComp() const {
   return _segComp;
}

float TCMuon::EmIso() const {
   return _emIso;
}

float TCMuon::HadIso() const {
   return _hadIso;
}

float TCMuon::TrkIso() const {
   return _trkIso;
}

int TCMuon::Ntracks() const {
   return _nTracks;
}

float TCMuon::PfSumPt(float coneSize) const {
  float sumPt = 0;
  if (fabs(coneSize - 0.3) < 0.01) sumPt = _pfIso_Pt03;
  if (fabs(coneSize - 0.4) < 0.01) sumPt = _pfIso_Pt04;
  if (fabs(coneSize - 0.5) < 0.01) sumPt = _pfIso_Pt05;
  return sumPt;
}

float TCMuon::PfENeutral(float coneSize) const {
  float neutral = 0;
  if (fabs(coneSize - 0.3) < 0.01) neutral = _pfIso_Neutral03;
  if (fabs(coneSize - 0.4) < 0.01) neutral = _pfIso_Neutral04;
  if (fabs(coneSize - 0.5) < 0.01) neutral = _pfIso_Neutral05;
  return neutral;
}

float TCMuon::PfEGamma(float coneSize) const {
  float gamma = 0;
  if (fabs(coneSize - 0.3) < 0.01) gamma = _pfIso_Gamma03;
  if (fabs(coneSize - 0.4) < 0.01) gamma = _pfIso_Gamma04;
  if (fabs(coneSize - 0.5) < 0.01) gamma = _pfIso_Gamma05;
  return gamma;
}

float TCMuon::PfRelIso(float coneSize) const {
  float relIso = 0;
  if (fabs(coneSize - 0.3) < 0.01)
    relIso = (_pfIso_Pt03 + _pfIso_Gamma03 + _pfIso_Neutral03) / _p4.Pt();
  if (fabs(coneSize - 0.4) < 0.01)
    relIso = (_pfIso_Pt04 + _pfIso_Gamma04 + _pfIso_Neutral04) / _p4.Pt();
  if (fabs(coneSize - 0.5) < 0.01)
    relIso = (_pfIso_Pt05 + _pfIso_Gamma05 + _pfIso_Neutral05) / _p4.Pt();
  return relIso;
}

/*
This is not working for some reason - check them later

float TCMuon::dxy(TVector3 *primVtx)
{
  //Calculating track dxy parameter wrt primary vertex
  //d0 = - dxy
  float vx = _vtx->X(), vy = _vtx->Y();
  float px = _p4->Px(), py = _p4->Py(), pt = _p4->Pt();
  float pvx = primVtx->X(), pvy = primVtx->Y();
  
  float ret =  (-(vx-pvx)*py + (vy-pvy)*px)/pt;
  return ret;
}

float TCMuon::dz(TVector3 *primVtx)
{
  //Calculating track dz parameter wrt primary vertex
  float vx = _vtx->X(), vy = _vtx->Y(), vz = _vtx->Z();
  float px = _p4->Px(), py = _p4->Py();
  float pz = _p4->Pz(), pt = _p4->Pt();
  float pvx = primVtx->X(), pvy = primVtx->Y(), pvz = primVtx->Z();
  
  float ret =  (vz-pvz)-((vx-pvx)*px +(vy-pvy)*py)/pt*(pz/pt);

  return ret;
}
*/

// "set" methods ---------------------------------------------

void TCMuon::SetP4(TLorentzVector p4) {
   _p4 = p4;
}

void TCMuon::SetP4(float px, float py, float pz, float e) {
   TLorentzVector p4(px, py, pz, e);
   _p4 = p4;
}

void TCMuon::SetVtx(float vx, float vy, float vz) {
   TVector3 v3(vx, vy, vz);
   _vtx = v3;
}

void TCMuon::SetPtError(float er){
   _ptError = er;
}
void TCMuon::SetCharge(int c){
   _charge = c;
}

void TCMuon::SetIsGLB(bool t){
   _isGLB = t;
}

void TCMuon::SetIsTRK(bool t){
   _isTRK = t;
}

void TCMuon::SetNumberOfValidMuonHits(int n) {
  _numberOfValidMuonHits = n;
}

void TCMuon::SetNumberOfValidPixelHits(int n) {
  _numberOfValidPixelHits = n;
}

void TCMuon::SetNumberOfValidTrackerHits(int n) {
  _numberOfValidTrackerHits = n;
}

void TCMuon::SetNumberOfLostPixelHits(int n) {
  _numberOfLostPixelHits = n;
}

void TCMuon::SetNumberOfLostTrackerHits(int n) {
  _numberOfLostTrackerHits = n;
}

void TCMuon::SetNormalizedChi2(float n) {
  _normalizedChi2 = n;
}

void TCMuon::SetNumberOfMatches(int n) {
  _numberOfMatches = n;
}

void TCMuon::SetCaloComp(float c){
   _caloComp = c;
}

void TCMuon::SetSegComp(float s){
   _segComp = s;
}

void TCMuon::SetEMIso(float e){
   _emIso = e;
}


void TCMuon::SetHADIso(float h){
   _hadIso = h;
}

void TCMuon::SetTRKIso(float t){
   _trkIso = t;
}

void TCMuon::SetNtracks(int n){
   _nTracks = n;
}
void TCMuon::SetPfSumPt(float coneSize, float f) {
  if(fabs(coneSize - 0.3) < 0.01) _pfIso_Pt03 = f;
  if(fabs(coneSize - 0.4) < 0.01) _pfIso_Pt04 = f;
  if(fabs(coneSize - 0.5) < 0.01) _pfIso_Pt05 = f;
}

void TCMuon::SetPfEGamma(float coneSize, float f) {
  if(fabs(coneSize - 0.3) < 0.01) _pfIso_Gamma03 = f;
  if(fabs(coneSize - 0.4) < 0.01) _pfIso_Gamma04 = f;
  if(fabs(coneSize - 0.5) < 0.01) _pfIso_Gamma05 = f;
}

void TCMuon::SetPfENeutral(float coneSize, float f) {
  if(fabs(coneSize - 0.3) < 0.01) _pfIso_Neutral03 = f;
  if(fabs(coneSize - 0.4) < 0.01) _pfIso_Neutral04 = f;
  if(fabs(coneSize - 0.5) < 0.01) _pfIso_Neutral05 = f;
}

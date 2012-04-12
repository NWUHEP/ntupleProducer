#include "Higgs/ntupleProducer/interface/TCMuon.h"
#include <iostream>

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

float TCMuon::EmIso03() const {
   return _emIso03;
}
float TCMuon::HadIso03() const {
   return _hadIso03;
}
float TCMuon::TrkIso03() const {
   return _trkIso03;
}


float TCMuon::EmIso() const {
   return _emIso03;
}
float TCMuon::HadIso() const {
   return _hadIso03;
}
float TCMuon::TrkIso() const {
   return _trkIso03;
}
float TCMuon::EmIso05() const {
   return _emIso05;
}
float TCMuon::HadIso05() const {
   return _hadIso05;
}
float TCMuon::TrkIso05() const {
   return _trkIso05;
}


int TCMuon::Ntracks() const {
   return _nTracks03;
}
int TCMuon::Ntracks03() const {
   return _nTracks03;
}
int TCMuon::Ntracks05() const {
   return _nTracks05;
}

float TCMuon::PfSumPt(int coneSize) const {
  float sumPt = 0;
  if (coneSize == 3) sumPt = _pfIso_Pt03;
  if (coneSize == 4) sumPt = _pfIso_Pt04;
  return sumPt;
}

float TCMuon::PfENeutral(int coneSize) const {
  float neutral = 0;
  if (coneSize == 3) neutral = _pfIso_Neutral03;
  if (coneSize == 4) neutral = _pfIso_Neutral04;
  return neutral;
}

float TCMuon::PfEGamma(int coneSize) const {
  float gamma = 0;
  if (coneSize == 3) gamma = _pfIso_Gamma03;
  if (coneSize == 4) gamma = _pfIso_Gamma04;
  return gamma;
}

float TCMuon::PfRelIso(int coneSize) const {
  float relIso = 0;
  if (coneSize == 3)
    relIso = (_pfIso_Pt03 + _pfIso_Gamma03 + _pfIso_Neutral03) / _p4.Pt();
  if (coneSize == 4)
    relIso = (_pfIso_Pt04 + _pfIso_Gamma04 + _pfIso_Neutral04) / _p4.Pt();
  return relIso;
}


float TCMuon::Dxy(TVector3 *primVtx) const {
  //Calculating track dxy parameter wrt primary vertex                                                                                                       
  //d0 = - dxy                                                                                                                                               
  float vx = _vtx.X(), vy = _vtx.Y();
  float px = _p4.Px(), py = _p4.Py(), pt = _p4.Pt();
  float pvx = primVtx->X(), pvy = primVtx->Y();
  float ret =  (-(vx-pvx)*py + (vy-pvy)*px)/pt;
  return ret;
}

float TCMuon::Dz(TVector3 *primVtx) const {
  //Calculating track dz parameter wrt primary vertex                                                                                                        
  float vx = _vtx.X(), vy = _vtx.Y(), vz = _vtx.Z();
  float px = _p4.Px(), py = _p4.Py();
  float pz = _p4.Pz(), pt = _p4.Pt();
  float pvx = primVtx->X(), pvy = primVtx->Y(), pvz = primVtx->Z();
  float ret =  (vz-pvz)-((vx-pvx)*px +(vy-pvy)*py)/pt*(pz/pt);
  return ret;
}

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


void TCMuon::SetEmIso03(float e){
   _emIso03 = e;
}
void TCMuon::SetHadIso03(float h){
   _hadIso03 = h;
}
void TCMuon::SetTrkIso03(float t){
   _trkIso03 = t;
}
void TCMuon::SetNtracks03(int n){
   _nTracks03 = n;
}


void TCMuon::SetEmIso05(float e){
   _emIso05 = e;
}
void TCMuon::SetHadIso05(float h){
   _hadIso05 = h;
}
void TCMuon::SetTrkIso05(float t){
   _trkIso05 = t;
}
void TCMuon::SetNtracks05(int n){
   _nTracks05 = n;
}

void TCMuon::SetPfSumPt(int coneSize, float f) {
  if(coneSize == 3) _pfIso_Pt03 = f;
  if(coneSize == 4) _pfIso_Pt04 = f;
}

void TCMuon::SetPfEGamma(int coneSize, float f) {
  if(coneSize == 3) _pfIso_Gamma03 = f;
  if(coneSize == 4) _pfIso_Gamma04 = f;
}

void TCMuon::SetPfENeutral(int coneSize, float f) {
  if(coneSize == 3) _pfIso_Neutral03 = f;
  if(coneSize == 4) _pfIso_Neutral04 = f;
}

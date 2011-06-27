#include "TCMuon.h"
#include<iostream>

TCMuon::TCMuon() {}

TCMuon::~TCMuon() {}


// "get" methods -------------------------------------

TLorentzVector TCMuon::p4() const {
   return _p4;
}

float TCMuon::pT() const {
   return _p4.Pt();
}

TVector3 TCMuon::Vtx() const {
   return _vtx;
}

//TVector3 TCMuon::AssocVtx() {
//   return _assocPV;
//}

float TCMuon::eta() const {
   return _eta;
}

float TCMuon::phi() const {
   return _phi;
}

int TCMuon::charge() const {
   return _charge;
}

bool TCMuon::isTRK() const {
   return _isTRK;
}

bool TCMuon::isGLB() const {
   return _isGLB;
}

float TCMuon::dxy() const {
   return _dxy;
}

int TCMuon::nPXLHits() const {
   return _nPXLHits;
}

int TCMuon::nTRKHits() const {
   return _nTRKHits;
}

int TCMuon::nMatchSeg() const {
   return _nMatchSeg;
}

int TCMuon::nValidMuHits() const {
   return _nValidMuHits;
}

float TCMuon::normChi2() const {
   return _normChi2;
}

float TCMuon::caloComp() const {
   return _caloComp;
}

float TCMuon::segComp() const {
   return _segComp;
}

float TCMuon::emIso() const {
   return _emIso;
}

float TCMuon::hadIso() const {
   return _hadIso;
}

float TCMuon::trkIso() const {
   return _trkIso;
}


// "set" methods ---------------------------------------------

void TCMuon::Setp4(TLorentzVector p4) {
   _p4 = p4;
}

void TCMuon::Setp4(float px, float py, float pz, float e) {
   TLorentzVector p4(px, py, pz, e);
   _p4 = p4;
}

void TCMuon::SetVtx(float vx, float vy, float vz) {
   TVector3 v3(vx, vy, vz);
   _vtx = v3;
}

//void TCMuon::SetAssocVtx(float vx, float vy, float vz) {
//   TVector3 v3(vx, vy, vz);
//   _assocPV = v3;
//}

void TCMuon::SetEta(float e){
   _eta = e;
}

void TCMuon::SetPhi(float p){
   _phi = p;
}

void TCMuon::SetCharge(int c){
   _charge = c;
}

void TCMuon::SetisGLB(bool t){
   _isGLB = t;
}

void TCMuon::SetisTRK(bool t){
   _isTRK = t;
}

void TCMuon::Setdxy(float d){
   _dxy = d;
}

void TCMuon::SetnPXLHits(int n){
   _nPXLHits = n;
}

void TCMuon::SetnMatchSeg(int n){
   _nMatchSeg = n;
}

void TCMuon::SetnTRKHits(int n){
   _nTRKHits = n;
}

void TCMuon::SetnValidMuHits(int n){
   _nValidMuHits = n;
}

void TCMuon::SetNormChi2(float c){
   _normChi2 = c;
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


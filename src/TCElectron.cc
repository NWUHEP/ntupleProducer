#include "TCElectron.h"
#include<iostream>

TCElectron::TCElectron() {
}

TCElectron::~TCElectron() {
}

// "get" methods -------------------------------------

TLorentzVector TCElectron::p4() const {
   return _p4;
}

float TCElectron::pT() const {
   return _p4.Pt();
}

TVector3 TCElectron::Vtx() const {
   return _vtx;
}

//TVector3 TCElectron::AssocVtx() {
//   return _assocPV;
//}

float TCElectron::eta() const {
   return _p4.Eta();
}

float TCElectron::phi() const {
   return _p4.Phi();
}

int TCElectron::charge() const {
   return _charge;
}

int TCElectron::idMap() const {
   return _idMap;
}

int TCElectron::conversionFlag() const {
   return _convFlag;
}

float TCElectron::dxy() const {
   return _dxy;
}

float TCElectron::normChi2() const {
   return _normChi2;
}

float TCElectron::emIso() const {
   return _emIso;
}

float TCElectron::hadIso() const {
   return _hadIso;
}

float TCElectron::trkIso() const {
   return _trkIso;
}

float TCElectron::HoverE() const {
	return _HoverE;
}

float TCElectron::dPhiSC() const {
	return _dPhiSC;
}

float TCElectron::dEtaSC() const {
	return _dEtaSC;
}

float TCElectron::Sig_IEtaIEta() const {
	return _sig_IEtaIEta;
}

float TCElectron::pfChargedHadronIso() const {
	return _pfChargedHadronIso;
}

float TCElectron::pfNeutralHadronIso() const {
	return _pfNeutralHadronIso;
}

float TCElectron::pfPhotonIso() const {
	return _pfPhotonIso;
}

float TCElectron::conversionDist() const {
	return _convDist;
}

float TCElectron::conversionDcot() const {
	return _convDcot;
}

float TCElectron::conversionRad() const {
	return _convRad;
}


// "set" methods ---------------------------------------------

void TCElectron::Setp4(TLorentzVector P4) {
   _p4 = P4;
}

void TCElectron::Setp4(float px, float py, float pz, float e) {
   TLorentzVector P4(px, py, pz, e);
   _p4 = P4;
}

void TCElectron::SetVtx(float vx, float vy, float vz) {
   TVector3 v3(vx, vy, vz);
   _vtx = v3;
}

//void TCElectron::SetAssocVtx(float vx, float vy, float vz) {
//   TVector3 v3(vx, vy, vz);
//   _assocPV = v3;
//}


void TCElectron::SetCharge(int c){
   _charge = c;
}

void TCElectron::SetIDMap(int i){
   _idMap = i;
}

void TCElectron::SetConversionFlag(int f){
   _convFlag = f;
}

void TCElectron::Setdxy(float d){
   _dxy = d;
}

void TCElectron::SetNormChi2(float c){
   _normChi2 = c;
}

void TCElectron::SetEMIso(float e){
   _emIso = e;
}

void TCElectron::SetHADIso(float h){
   _hadIso = h;
}

void TCElectron::SetTRKIso(float t){
   _trkIso = t;
}

void TCElectron::SetHoverE(float h){
	_HoverE = h;
}

void TCElectron::SetdPhiSC(float d){
	_dPhiSC = d;
}

void TCElectron::SetdEtaSC(float d){
	_dEtaSC = d;
}

void TCElectron::SetSig_IEtaIEta(float s){
	_sig_IEtaIEta = s;
}

void TCElectron::SetPFChargedHadronIso(float c){
	_pfChargedHadronIso = c;
}

void TCElectron::SetPFNeutralHadronIso(float n){
	_pfNeutralHadronIso = n;
}

void TCElectron::SetPFPhotonIso(float g){
	_pfPhotonIso = g;
}

void TCElectron::SetConversionDist(float d) {
	_convDist = d;
}

void TCElectron::SetConversionDcot(float d) {
	_convDcot = d;
}

void TCElectron::SetConversionRad(float r) {
	_convRad = r;
}


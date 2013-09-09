#include "../interface/TCMuon.h"
#include "../interface/TCMuonLinkDef.h"
#include <iostream>

TCMuon::TCMuon() {
}
TCMuon::~TCMuon() {
}


// "get" methods -------------------------------------


float TCMuon::PtError() const {
   return _ptError;
}

bool TCMuon::IsGLB() const {
   return _isGLB;
}
bool TCMuon::IsTRK() const {
   return _isTRK;
}

bool TCMuon::IsPF() const {
   return _isPF;
}

int TCMuon::NumberOfMatchedStations() const {
  return _numberOfMatchedStations;
}

int TCMuon::TrackLayersWithMeasurement() const {
  return _trackLayersWithMeasurement;
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

// "set" methods ---------------------------------------------
void TCMuon::SetNumberOfMatchedStations(int n){
  _numberOfMatchedStations = n;
}

void TCMuon::SetTrackLayersWithMeasurement(int n) {
  _trackLayersWithMeasurement = n;
}

void TCMuon::SetPtError(float er){
   _ptError = er;
}

void TCMuon::SetIsGLB(bool t){
   _isGLB = t;
}

void TCMuon::SetIsTRK(bool t){
   _isTRK = t;
}
void TCMuon::SetIsPF(bool t){
   _isPF = t;
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


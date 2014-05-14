#include "../interface/TCMuon.h"
#include "TCMuonLinkDef.h"

TCMuon::TCMuon():
_ptError(-99), 
_isTRK(0),
_isGLB(0),
_isSoft(0),
_isTight(0),
_isGood(0),
_isGoodLoose(0),
_isArbitrated(0),
_isTrkArbitrated(0),
_caloComp(-99),
_segComp(-99),
_numberOfMatches(-99),
_numberOfMatchedStations(-99),
_numberOfValidPixelHits(-99),
_numberOfValidTrackerHits(-99),
_numberOfLostPixelHits(-99),
_numberOfLostTrackerHits(-99),
_numberOfValidMuonHits(-99),
_trackLayersWithMeasurement(-99),
_pixelLayersWithMeasurement(-99),
_normalizedChi2(-99),
_normalizedChi2_tracker(-99),
_pfIsoPU(-99),
_pfIsoChargedHad(-99),
_pfIsoChargedPart(-99),
_pfIsoNeutral(-99),
_pfIsoPhoton(-99)
{
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

bool TCMuon::IsSoft() const {
  return _isSoft;
}
bool TCMuon::IsTight() const {
  return _isTight;
}

bool TCMuon::IsGood() const {
  return _isGood;
}
bool TCMuon::IsGoodLoose() const {
  return _isGoodLoose;
}

bool TCMuon::IsArbitrated() const {
  return _isArbitrated;
}
bool TCMuon::IsTrkArbitrated() const {
  return _isTrkArbitrated;
}

int TCMuon::NumberOfMatchedStations() const {
  return _numberOfMatchedStations;
}

int TCMuon::TrackLayersWithMeasurement() const {
  return _trackLayersWithMeasurement;
}
int TCMuon::PixelLayersWithMeasurement() const {
  return _pixelLayersWithMeasurement;
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
float TCMuon::NormalizedChi2_tracker() const {
  return _normalizedChi2_tracker;
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

float TCMuon::PfIsoPU() const {
  return _pfIsoPU;
}
float TCMuon::PfIsoCharged() const {
  return _pfIsoChargedHad;
}
float TCMuon::PfIsoChargedHad() const {
  return _pfIsoChargedHad;
}
float TCMuon::PfIsoChargedPart() const {
  return _pfIsoChargedPart;
}
float TCMuon::PfIsoNeutral() const {
  return _pfIsoNeutral;
}
float TCMuon::PfIsoPhoton() const {
  return _pfIsoPhoton;
}


// "set" methods ---------------------------------------------
void TCMuon::SetNumberOfMatchedStations(int n){
  _numberOfMatchedStations = n;
}

void TCMuon::SetTrackLayersWithMeasurement(int n) {
  _trackLayersWithMeasurement = n;
}
void TCMuon::SetPixelLayersWithMeasurement(int n) {
  _pixelLayersWithMeasurement = n;
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

void TCMuon::SetIsSoft(bool t){
  _isSoft = t;
}
void TCMuon::SetIsTight(bool t){
  _isTight = t;
}
void TCMuon::SetIsGood(bool g){
  _isGood = g;
}
void TCMuon::SetIsGoodLoose(bool g){
  _isGoodLoose = g;
}
void TCMuon::SetIsArbitrated(bool g){
  _isArbitrated = g;
}
void TCMuon::SetIsTrkArbitrated(bool g){
  _isTrkArbitrated= g;
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
void TCMuon::SetNormalizedChi2_tracker(float n) {
  _normalizedChi2_tracker = n;
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

void TCMuon::SetPfIsoPU(float f) {
  _pfIsoPU = f;
}
void TCMuon::SetPfIsoChargedHad(float f) {
  _pfIsoChargedHad = f;
}
void TCMuon::SetPfIsoChargedPart(float f) {
  _pfIsoChargedPart = f;
}
void TCMuon::SetPfIsoNeutral(float f) {
  _pfIsoNeutral = f;
}
void TCMuon::SetPfIsoPhoton(float f) {
  _pfIsoPhoton = f;
}


ostream& TCMuon::TCprint(ostream& os) const {
 return TCPhysObject::TCprint(os) << 
   " IsPF: "<< IsPF() << " IsGLB: " << IsGLB() << " IsTRK: " << IsTRK() << " IsSoft: "<< IsSoft() << " IsTight: "<< IsTight() <<
   " IsGood: " << IsGood() << " IsGoodLoose: " << IsGoodLoose() << " CaloComp: " << CaloComp() << " SegComp: " << SegComp() <<
   " NumberOfValidPixelHits: " << NumberOfValidPixelHits() << " NumberOfValidTrackerHits: " << NumberOfValidTrackerHits() <<
   " NumberOfValidMuonHits: " << NumberOfValidMuonHits() << " NumberOfLostPixelHits: " << NumberOfLostPixelHits() << 
   " NumberOfLostTrackerHits: " << NumberOfLostTrackerHits() << " NumberOfMatches: " << NumberOfMatches() << 
   " NumberOfMatchedStations: " << NumberOfMatchedStations() << " TrackLayersWithMeasurement: " << TrackLayersWithMeasurement() <<
   " PixelLayersWithMeasurement: " << PixelLayersWithMeasurement() << " NormalizedChi2: " << NormalizedChi2() << 
   " NormalizedChi2_tracker: " << NormalizedChi2_tracker() << " PfIsoPU: " << PfIsoPU() << " PfIsoChargedPart: " << PfIsoChargedPart() <<
   " PfIsoChargedHad: " << PfIsoChargedHad() << " PfIsoCharged: " << PfIsoCharged() << " PfIsoNeutral: " << PfIsoNeutral() << 
   " PfIsoPhoton: " << PfIsoPhoton();
}


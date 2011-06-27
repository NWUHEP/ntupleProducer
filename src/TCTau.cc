#include "TCTau.h"
#include <iostream>

TCTau::TCTau() {
}

TCTau::~TCTau() {
}

TVector3 TCTau::PositionFromTrack() const {
	   return _positionFromLeadTrack;
}

TVector3 TCTau::PositionFromTau() const {
	   return _positionFromTauObject;
}

TLorentzVector TCTau::P4() const {
	   return _p4;
}

TVector2 TCTau::P2() const {
	   TVector2 v2(_p4.Px(), _p4.Py());
		   return v2;
}

float TCTau::Et() const {
	   return _p4.Et();
}

float TCTau::Pt() const {
	   return _p4.Pt();
}

int TCTau::GetCharge() {
	   return pftau_charge;
}

int TCTau::GetNTracks() {
	   return ntracks;
}

int TCTau::GetNConst() {
	   return nconst;
}

std::vector<TLorentzVector> TCTau::GetChargedTracks() {
	   return pftau_ch;
}

void TCTau::AddChargedTrack(TLorentzVector v){
	   pftau_ch.push_back(v);
}

std::vector<TLorentzVector> TCTau::GetH0Tracks() {
	   return pftau_h0;
}

void TCTau::AddH0Track(TLorentzVector v){
	   pftau_h0.push_back(v);
}

std::vector<TLorentzVector> TCTau::GetEMTracks() {
	   return pftau_em;
}

void TCTau::AddEMTrack(TLorentzVector v){
	   pftau_em.push_back(v);
}

void TCTau::SetCharge(int c) {
	   pftau_charge = c;
}

void TCTau::SetNTracks(int n) {
	   ntracks = n;
}

void TCTau::SetNConst(int n) {
	   nconst = n;
}

void TCTau::SetleadPFChargedPt(float p) {
	   leadPFCharged = p;
}

void TCTau::SetleadPFNeutralPt(float p) {
	   leadPFNeutral = p;
}

float TCTau::leadPFChargedPt() {
	   return leadPFCharged;
}

float TCTau::leadPFNeutralPt() {
	   return leadPFNeutral;
}

float TCTau::leadPFAnyPt() {
	   if(leadPFNeutral > leadPFCharged) {
			      return leadPFNeutral;   
					   }else{
							      return leadPFCharged;
									   }
}

void TCTau::SetDecayMode(int d) {
	   DecayMode = d;
}

int TCTau::GetDecayMode() {
	   return DecayMode;
}

void TCTau::SetPositionFromTau(float x, float y, float z) {
	   TVector3 p(x, y, z);
		   _positionFromTauObject = p;
}

void TCTau::SetPositionFromTrack(float x, float y, float z) {
	   TVector3 p(x, y, z);
		   _positionFromLeadTrack = p;
}

void TCTau::SetP4(TLorentzVector p4) {
	   _p4 = p4;
}

void TCTau::SetP4(float px, float py, float pz, float e) {
	   TLorentzVector p4(px, py, pz, e);
		   _p4 = p4;
}

void TCTau::SetAlgorithm(int i) {
	   TauAlgorithm = i;
}

int TCTau::GetAlgorithm() const{
	  return TauAlgorithm;
}

void TCTau::SethpsPFTauDiscriminationAgainstElectronLoose(float f) {
	   hpsPFTauDiscriminationAgainstElectronLoose = f;
}

void TCTau::SethpsPFTauDiscriminationAgainstMuonLoose(float f) {
	   hpsPFTauDiscriminationAgainstMuonLoose = f;
}

void TCTau::SethpsPFTauDiscriminationAgainstElectronMedium(float f) {
	   hpsPFTauDiscriminationAgainstElectronMedium = f;
}
   
void TCTau::SethpsPFTauDiscriminationAgainstMuonMedium(float f) {
	   hpsPFTauDiscriminationAgainstMuonMedium = f;
}

void TCTau::SethpsPFTauDiscriminationAgainstElectronTight(float f) {
	   hpsPFTauDiscriminationAgainstElectronTight = f;
}
   
void TCTau::SethpsPFTauDiscriminationAgainstMuonTight(float f) {
	   hpsPFTauDiscriminationAgainstMuonTight = f;
}

void TCTau::SethpsPFTauDiscriminationByDecayModeFinding(float f) {
	   hpsPFTauDiscriminationByDecayModeFinding = f;
}

void TCTau::SethpsPFTauDiscriminationByLooseIsolation(float f) {
	   hpsPFTauDiscriminationByLooseIsolation = f;
}

void TCTau::SethpsPFTauDiscriminationByMediumIsolation(float f) {
	   hpsPFTauDiscriminationByMediumIsolation = f;
}

void TCTau::SethpsPFTauDiscriminationByTightIsolation(float f) {
	   hpsPFTauDiscriminationByTightIsolation = f;
}

void TCTau::SetshrinkingConePFTauDecayModeIndexProducer(float f) {
	   shrinkingConePFTauDecayModeIndexProducer = f;
}

void TCTau::SetshrinkingConePFTauDiscriminationAgainstElectron(float f) {
	   shrinkingConePFTauDiscriminationAgainstElectron = f;
}

void TCTau::SetshrinkingConePFTauDiscriminationAgainstMuon(float f) {
	   shrinkingConePFTauDiscriminationAgainstMuon = f;
}

void TCTau::SetshrinkingConePFTauDiscriminationByECALIsolation(float f) {
	   shrinkingConePFTauDiscriminationByECALIsolation = f;
}

void TCTau::SetshrinkingConePFTauDiscriminationByECALIsolationUsingLeadingPion(float f) {
	   shrinkingConePFTauDiscriminationByECALIsolationUsingLeadingPion = f;
}

void TCTau::SetshrinkingConePFTauDiscriminationByIsolation(float f) {
	   shrinkingConePFTauDiscriminationByIsolation = f;
}

void TCTau::SetshrinkingConePFTauDiscriminationByIsolationUsingLeadingPion(float f) {
	   shrinkingConePFTauDiscriminationByIsolationUsingLeadingPion = f;
}

void TCTau::SetshrinkingConePFTauDiscriminationByLeadingPionPtCut(float f) {
	   shrinkingConePFTauDiscriminationByLeadingPionPtCut = f;
}

void TCTau::SetshrinkingConePFTauDiscriminationByLeadingTrackFinding(float f) {
	   shrinkingConePFTauDiscriminationByLeadingTrackFinding = f;
}

void TCTau::SetshrinkingConePFTauDiscriminationByLeadingTrackPtCut(float f) {
	   shrinkingConePFTauDiscriminationByLeadingTrackPtCut = f;
}

void TCTau::SetshrinkingConePFTauDiscriminationByTaNC(float f) {
	   shrinkingConePFTauDiscriminationByTaNC = f;
}

void TCTau::SetshrinkingConePFTauDiscriminationByTaNCfrHalfPercent(float f) {
	   shrinkingConePFTauDiscriminationByTaNCfrHalfPercent = f;
}

void TCTau::SetshrinkingConePFTauDiscriminationByTaNCfrOnePercent(float f) {
	   shrinkingConePFTauDiscriminationByTaNCfrOnePercent = f;
}

void TCTau::SetshrinkingConePFTauDiscriminationByTaNCfrQuarterPercent(float f) {
	   shrinkingConePFTauDiscriminationByTaNCfrQuarterPercent = f;
}

void TCTau::SetshrinkingConePFTauDiscriminationByTaNCfrTenthPercent(float f) {
	   shrinkingConePFTauDiscriminationByTaNCfrTenthPercent = f;
}

void TCTau::SetshrinkingConePFTauDiscriminationByTrackIsolation(float f) {
	   shrinkingConePFTauDiscriminationByTrackIsolation = f;
}

void TCTau::SetshrinkingConePFTauDiscriminationByTrackIsolationUsingLeadingPion(float f) {
	   shrinkingConePFTauDiscriminationByTrackIsolationUsingLeadingPion = f;
}

float TCTau::GethpsPFTauDiscriminationAgainstElectronLoose() {
	   return hpsPFTauDiscriminationAgainstElectronLoose;
}

float TCTau::GethpsPFTauDiscriminationAgainstMuonLoose() {
	   return hpsPFTauDiscriminationAgainstMuonLoose;
}
   
float TCTau::GethpsPFTauDiscriminationAgainstElectronMedium() {
	   return hpsPFTauDiscriminationAgainstElectronMedium;
}
   
float TCTau::GethpsPFTauDiscriminationAgainstMuonMedium() {
	   return hpsPFTauDiscriminationAgainstMuonMedium;
}
   
float TCTau::GethpsPFTauDiscriminationAgainstElectronTight() {
	   return hpsPFTauDiscriminationAgainstElectronTight;
}
   
float TCTau::GethpsPFTauDiscriminationAgainstMuonTight() {
	   return hpsPFTauDiscriminationAgainstMuonTight;
}

float TCTau::GethpsPFTauDiscriminationByDecayModeFinding() {
	   return hpsPFTauDiscriminationByDecayModeFinding;
}

float TCTau::GethpsPFTauDiscriminationByLooseIsolation() {
	   return hpsPFTauDiscriminationByLooseIsolation;
}

float TCTau::GethpsPFTauDiscriminationByMediumIsolation() {
	   return hpsPFTauDiscriminationByMediumIsolation;
}

float TCTau::GethpsPFTauDiscriminationByTightIsolation() {
	   return hpsPFTauDiscriminationByTightIsolation;
}

float TCTau::GetshrinkingConePFTauDecayModeIndexProducer() {
	   return shrinkingConePFTauDecayModeIndexProducer;
}

float TCTau::GetshrinkingConePFTauDiscriminationAgainstElectron() {
	   return shrinkingConePFTauDiscriminationAgainstElectron;
}

float TCTau::GetshrinkingConePFTauDiscriminationAgainstMuon() {
	   return shrinkingConePFTauDiscriminationAgainstMuon;
}

float TCTau::GetshrinkingConePFTauDiscriminationByECALIsolation() {
	   return shrinkingConePFTauDiscriminationByECALIsolation;
}

float TCTau::GetshrinkingConePFTauDiscriminationByECALIsolationUsingLeadingPion() {
	   return shrinkingConePFTauDiscriminationByECALIsolationUsingLeadingPion;
}

float TCTau::GetshrinkingConePFTauDiscriminationByIsolation() {
	   return shrinkingConePFTauDiscriminationByIsolation;
}

float TCTau::GetshrinkingConePFTauDiscriminationByIsolationUsingLeadingPion() {
	   return shrinkingConePFTauDiscriminationByIsolationUsingLeadingPion;
}

float TCTau::GetshrinkingConePFTauDiscriminationByLeadingPionPtCut() {
	   return shrinkingConePFTauDiscriminationByLeadingPionPtCut;
}

float TCTau::GetshrinkingConePFTauDiscriminationByLeadingTrackFinding() {
	   return shrinkingConePFTauDiscriminationByLeadingTrackFinding;
}

float TCTau::GetshrinkingConePFTauDiscriminationByLeadingTrackPtCut() {
	   return shrinkingConePFTauDiscriminationByLeadingTrackPtCut;
}

float TCTau::GetshrinkingConePFTauDiscriminationByTaNC() {
	   return shrinkingConePFTauDiscriminationByTaNC;
}

float TCTau::GetshrinkingConePFTauDiscriminationByTaNCfrHalfPercent() {
	   return shrinkingConePFTauDiscriminationByTaNCfrHalfPercent;
}

float TCTau::GetshrinkingConePFTauDiscriminationByTaNCfrOnePercent() {
	   return shrinkingConePFTauDiscriminationByTaNCfrOnePercent;
}

float TCTau::GetshrinkingConePFTauDiscriminationByTaNCfrQuarterPercent() {
	   return shrinkingConePFTauDiscriminationByTaNCfrQuarterPercent;
}

float TCTau::GetshrinkingConePFTauDiscriminationByTaNCfrTenthPercent() {
	   return shrinkingConePFTauDiscriminationByTaNCfrTenthPercent;
}

float TCTau::GetshrinkingConePFTauDiscriminationByTrackIsolation() {
	   return shrinkingConePFTauDiscriminationByTrackIsolation;
}

float TCTau::GetshrinkingConePFTauDiscriminationByTrackIsolationUsingLeadingPion() {
	   return shrinkingConePFTauDiscriminationByTrackIsolationUsingLeadingPion;
}


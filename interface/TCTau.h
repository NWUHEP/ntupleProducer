#ifndef _TCTau_H
#define _TCTau_H

#include "TObject.h"
#include "TLorentzVector.h"
#include <vector>

class TCTau : public TObject {
	private:
		TVector3 _positionFromTauObject;
		TVector3 _positionFromLeadTrack;
		TLorentzVector _p4;

		std::vector<TLorentzVector> pftau_ch;
		std::vector<TLorentzVector> pftau_h0;
		std::vector<TLorentzVector> pftau_em;
		int pftau_charge;
		int ntracks;
		int nconst;

		float hpsPFTauDiscriminationAgainstElectronLoose;        
		float hpsPFTauDiscriminationAgainstMuonLoose;        
		float hpsPFTauDiscriminationAgainstElectronMedium;        
		float hpsPFTauDiscriminationAgainstMuonMedium;        
		float hpsPFTauDiscriminationAgainstElectronTight;        
		float hpsPFTauDiscriminationAgainstMuonTight;        
		float hpsPFTauDiscriminationByDecayModeFinding;        
		float hpsPFTauDiscriminationByLooseIsolation;        
		float hpsPFTauDiscriminationByMediumIsolation;        
		float hpsPFTauDiscriminationByTightIsolation;        
		float shrinkingConePFTauDecayModeIndexProducer;        
		float shrinkingConePFTauDiscriminationAgainstElectron;        
		float shrinkingConePFTauDiscriminationAgainstMuon;        
		float shrinkingConePFTauDiscriminationByECALIsolation;        
		float shrinkingConePFTauDiscriminationByECALIsolationUsingLeadingPion;        
		float shrinkingConePFTauDiscriminationByIsolation;        
		float shrinkingConePFTauDiscriminationByIsolationUsingLeadingPion;        
		float shrinkingConePFTauDiscriminationByLeadingPionPtCut;        
		float shrinkingConePFTauDiscriminationByLeadingTrackFinding;        
		float shrinkingConePFTauDiscriminationByLeadingTrackPtCut;        
		float shrinkingConePFTauDiscriminationByTaNC;        
		float shrinkingConePFTauDiscriminationByTaNCfrHalfPercent;        
		float shrinkingConePFTauDiscriminationByTaNCfrOnePercent;        
		float shrinkingConePFTauDiscriminationByTaNCfrQuarterPercent;        
		float shrinkingConePFTauDiscriminationByTaNCfrTenthPercent;        
		float shrinkingConePFTauDiscriminationByTrackIsolation;        
		float shrinkingConePFTauDiscriminationByTrackIsolationUsingLeadingPion;        

		int DecayMode;
		float leadPFCharged;
		float leadPFNeutral;

		int TauAlgorithm;

	public:
		TCTau();
		virtual ~TCTau();

		TVector3 PositionFromTrack() const;
		TVector3 PositionFromTau() const;
		TLorentzVector P4() const;
		TVector2 P2() const;
		float Et() const;
		float Pt() const;
		int GetAlgorithm() const;

		float GethpsPFTauDiscriminationAgainstElectronLoose();
		float GethpsPFTauDiscriminationAgainstMuonLoose();
		float GethpsPFTauDiscriminationAgainstElectronMedium();
		float GethpsPFTauDiscriminationAgainstMuonMedium();
		float GethpsPFTauDiscriminationAgainstElectronTight();
		float GethpsPFTauDiscriminationAgainstMuonTight();
		float GethpsPFTauDiscriminationByDecayModeFinding();
		float GethpsPFTauDiscriminationByLooseIsolation();
		float GethpsPFTauDiscriminationByMediumIsolation();
		float GethpsPFTauDiscriminationByTightIsolation();
		float GetshrinkingConePFTauDecayModeIndexProducer();
		float GetshrinkingConePFTauDiscriminationAgainstElectron();
		float GetshrinkingConePFTauDiscriminationAgainstMuon();
		float GetshrinkingConePFTauDiscriminationByECALIsolation();
		float GetshrinkingConePFTauDiscriminationByECALIsolationUsingLeadingPion();
		float GetshrinkingConePFTauDiscriminationByIsolation();
		float GetshrinkingConePFTauDiscriminationByIsolationUsingLeadingPion();
		float GetshrinkingConePFTauDiscriminationByLeadingPionPtCut();
		float GetshrinkingConePFTauDiscriminationByLeadingTrackFinding();
		float GetshrinkingConePFTauDiscriminationByLeadingTrackPtCut();
		float GetshrinkingConePFTauDiscriminationByTaNC();
		float GetshrinkingConePFTauDiscriminationByTaNCfrHalfPercent();
		float GetshrinkingConePFTauDiscriminationByTaNCfrOnePercent();
		float GetshrinkingConePFTauDiscriminationByTaNCfrQuarterPercent();
		float GetshrinkingConePFTauDiscriminationByTaNCfrTenthPercent();
		float GetshrinkingConePFTauDiscriminationByTrackIsolation();
		float GetshrinkingConePFTauDiscriminationByTrackIsolationUsingLeadingPion();

		std::vector<TLorentzVector> GetChargedTracks();
		std::vector<TLorentzVector> GetH0Tracks();
		std::vector<TLorentzVector> GetEMTracks();
		int GetDecayMode();
		float leadPFChargedPt();
		float leadPFNeutralPt();
		float leadPFAnyPt();

		int GetCharge();
		int GetNTracks();
		int GetNConst();

		void AddChargedTrack(TLorentzVector v);
		void AddEMTrack(TLorentzVector v);
		void AddH0Track(TLorentzVector v);

		void SethpsPFTauDiscriminationAgainstElectronLoose(float f);
		void SethpsPFTauDiscriminationAgainstMuonLoose(float f);
		void SethpsPFTauDiscriminationAgainstElectronMedium(float f);
		void SethpsPFTauDiscriminationAgainstMuonMedium(float f);
		void SethpsPFTauDiscriminationAgainstElectronTight(float f);
		void SethpsPFTauDiscriminationAgainstMuonTight(float f);
		void SethpsPFTauDiscriminationByDecayModeFinding(float f);
		void SethpsPFTauDiscriminationByLooseIsolation(float f);
		void SethpsPFTauDiscriminationByMediumIsolation(float f);
		void SethpsPFTauDiscriminationByTightIsolation(float f);
		void SetshrinkingConePFTauDecayModeIndexProducer(float f);
		void SetshrinkingConePFTauDiscriminationAgainstElectron(float f);
		void SetshrinkingConePFTauDiscriminationAgainstMuon(float f);
		void SetshrinkingConePFTauDiscriminationByECALIsolation(float f);
		void SetshrinkingConePFTauDiscriminationByECALIsolationUsingLeadingPion(float f);
		void SetshrinkingConePFTauDiscriminationByIsolation(float f);
		void SetshrinkingConePFTauDiscriminationByIsolationUsingLeadingPion(float f);
		void SetshrinkingConePFTauDiscriminationByLeadingPionPtCut(float f);
		void SetshrinkingConePFTauDiscriminationByLeadingTrackFinding(float f);
		void SetshrinkingConePFTauDiscriminationByLeadingTrackPtCut(float f);
		void SetshrinkingConePFTauDiscriminationByTaNC(float f);
		void SetshrinkingConePFTauDiscriminationByTaNCfrHalfPercent(float f);
		void SetshrinkingConePFTauDiscriminationByTaNCfrOnePercent(float f);
		void SetshrinkingConePFTauDiscriminationByTaNCfrQuarterPercent(float f);
		void SetshrinkingConePFTauDiscriminationByTaNCfrTenthPercent(float f);
		void SetshrinkingConePFTauDiscriminationByTrackIsolation(float f);
		void SetshrinkingConePFTauDiscriminationByTrackIsolationUsingLeadingPion(float f);

		void SetleadPFChargedPt(float p);
		void SetleadPFNeutralPt(float p);
		void SetDecayMode(int d);
		void SetCharge(int c);
		void SetNTracks(int n);
		void SetNConst(int n);

		void SetPositionFromTau(float x, float y, float z);
		void SetPositionFromTrack(float x, float y, float z);
		void SetP4(TLorentzVector p4);
		void SetP4(float px, float py, float pz, float e);
		void SetAlgorithm(int i);

		ClassDef(TCTau, 4);

};

#endif

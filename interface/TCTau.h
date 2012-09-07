#ifndef _TCTau_H
#define _TCTau_H

#include "TObject.h"
#include "TLorentzVector.h"
#include "TCPhysObject.h"
#include <vector>

class TCTau : public TCPhysObject {
	private:
		TVector3 _positionFromTauObject;
		TVector3 _positionFromLeadTrack;

		int nChHad;
		int nGamma;
		int nNeutrHad;
		int decayMode;

		TLorentzVector leadChHadP4;
		TLorentzVector leadNeutrP4;

		float isoGammaEtSum;
		float isoChHadPtSum;


		// Discriminators - names as used in CMSSW
		// For the most part the discriminators are binary

		// HPS discriminators ---------
		//
		// against e,mu misidentification
		float hpsPFTauDiscriminationAgainstElectronLoose;        
		float hpsPFTauDiscriminationAgainstMuonLoose;        
		float hpsPFTauDiscriminationAgainstElectronMedium;        
		float hpsPFTauDiscriminationAgainstMuonMedium;        
		float hpsPFTauDiscriminationAgainstElectronTight;        
		float hpsPFTauDiscriminationAgainstMuonTight;        
		
		// by decay mode 
		float hpsPFTauDiscriminationByDecayModeFinding;

		// by charged isolation
		float hpsPFTauDiscriminationByVLooseIsolation;
		float hpsPFTauDiscriminationByLooseIsolation;        
		float hpsPFTauDiscriminationByMediumIsolation;        
		float hpsPFTauDiscriminationByTightIsolation;
		
		// with iso corrections
		float hpsPFTauDiscriminationByVLooseIsolationDBSumPtCorr;
		float hpsPFTauDiscriminationByLooseIsolationDBSumPtCorr;
		float hpsPFTauDiscriminationByMediumIsolationDBSumPtCorr;
		float hpsPFTauDiscriminationByTightIsolationDBSumPtCorr;

		// by combined isolation with iso corrections
		float hpsPFTauDiscriminationByVLooseCombinedIsolationDBSumPtCorr;
		float hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr;
		float hpsPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr;
		float hpsPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr;


	public:
		TCTau();
		virtual ~TCTau();

		TVector3 PositionFromTrack() const;
		TVector3 PositionFromTau() const;

		int DecayMode() const { return decayMode; }
		int NChHad() const  {return nChHad; }
		int NGamma() const {return nGamma; }
		int NNeutrHad() const {return nNeutrHad; }

		TLorentzVector LeadChHadP4() const { return leadChHadP4; }
		TLorentzVector LeadNeutrP4() const { return leadNeutrP4; }

		float IsoGammaEtSum() { return isoGammaEtSum; }
		float IsoChHadPtSum() { return isoChHadPtSum; }

		void SetDecayMode(int d) {decayMode = d; } 
		void SetNChHad(int n) { nChHad = n; }
		void SetNGamma(int n) { nGamma = n; }
		void SetNNeutrHad(int n) { nNeutrHad = n; }

		void SetPositionFromTau(float x, float y, float z);
		void SetPositionFromTrack(float x, float y, float z);

		void SetLeadChHadP4(TLorentzVector p4) { leadChHadP4 = p4; }
		void SetLeadNeutrP4(TLorentzVector p4) { leadNeutrP4 = p4; }

		void SetLeadChHadP4(float x, float y, float z, float e) { leadChHadP4.SetPxPyPzE(x,y,z,e); }
		void SetLeadNeutrP4(float x, float y, float z, float e) { leadNeutrP4.SetPxPyPzE(x,y,z,e); }

		void SetIsoGammaEtSum(float i) { isoGammaEtSum = i; }
		void SetIsoChHadPtSum(float i) { isoChHadPtSum = i; }


		// Discriminators: for now only for HPS

		// HPS
		float GetHpsPFTauDiscriminationAgainstElectronLoose() {
		  return hpsPFTauDiscriminationAgainstElectronLoose;
		}
		float GetHpsPFTauDiscriminationAgainstMuonLoose() {
		  return hpsPFTauDiscriminationAgainstMuonLoose;
		}
		float GetHpsPFTauDiscriminationAgainstElectronMedium() {
		  return hpsPFTauDiscriminationAgainstElectronMedium;
		}
		float GetHpsPFTauDiscriminationAgainstMuonMedium() {
		  return hpsPFTauDiscriminationAgainstMuonMedium;
		}
		float GetHpsPFTauDiscriminationAgainstElectronTight() {
		  return hpsPFTauDiscriminationAgainstElectronTight;
		}
		float GetHpsPFTauDiscriminationAgainstMuonTight() {
		  return hpsPFTauDiscriminationAgainstMuonTight;
		}
		float GetHpsPFTauDiscriminationByDecayModeFinding() {
		  return hpsPFTauDiscriminationByDecayModeFinding;
		}
		float GetHpsPFTauDiscriminationByVLooseIsolation() {
		  return hpsPFTauDiscriminationByVLooseIsolation;
		}
		float GetHpsPFTauDiscriminationByLooseIsolation() {
		  return hpsPFTauDiscriminationByLooseIsolation;
		}
		float GetHpsPFTauDiscriminationByMediumIsolation() {
		  return hpsPFTauDiscriminationByMediumIsolation;
		}
		float GetHpsPFTauDiscriminationByTightIsolation() {
		  return hpsPFTauDiscriminationByTightIsolation;
		}
		


		float GetHpsPFTauDiscriminationByVLooseIsolationDBSumPtCorr() {
		  return hpsPFTauDiscriminationByVLooseIsolationDBSumPtCorr;
		}
		float GetHpsPFTauDiscriminationByLooseIsolationDBSumPtCorr() {
		  return hpsPFTauDiscriminationByLooseIsolationDBSumPtCorr;
		}
		float GetHpsPFTauDiscriminationByMediumIsolationDBSumPtCorr() {
		  return hpsPFTauDiscriminationByMediumIsolationDBSumPtCorr;
		}
		float GetHpsPFTauDiscriminationByTightIsolationDBSumPtCorr() {
		  return hpsPFTauDiscriminationByTightIsolationDBSumPtCorr;
		}


		float GetHpsPFTauDiscriminationByVLooseCombinedIsolationDBSumPtCorr() {
		  return hpsPFTauDiscriminationByVLooseCombinedIsolationDBSumPtCorr;
		}
		float GetHpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr() {
		  return hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr;
		}
		float GetHpsPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr() {
		  return hpsPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr;
		}
		float GetHpsPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr() {
		  return hpsPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr;
		}



		// HPS discriminator setters
		void SetHpsPFTauDiscriminationAgainstElectronLoose(float f) {
		  hpsPFTauDiscriminationAgainstElectronLoose = f;
		}
		void SetHpsPFTauDiscriminationAgainstMuonLoose(float f) {
		  hpsPFTauDiscriminationAgainstMuonLoose = f;
		}
		void SetHpsPFTauDiscriminationAgainstElectronMedium(float f) {
		  hpsPFTauDiscriminationAgainstElectronMedium = f;
		}
		void SetHpsPFTauDiscriminationAgainstMuonMedium(float f) {
		  hpsPFTauDiscriminationAgainstMuonMedium = f;
		}
		void SetHpsPFTauDiscriminationAgainstElectronTight(float f) {
		  hpsPFTauDiscriminationAgainstElectronTight = f;
		}
		void SetHpsPFTauDiscriminationAgainstMuonTight(float f) {
		  hpsPFTauDiscriminationAgainstMuonTight = f;
		}
		void SetHpsPFTauDiscriminationByDecayModeFinding(float f) {
		  hpsPFTauDiscriminationByDecayModeFinding = f;
		}
		void SetHpsPFTauDiscriminationByVLooseIsolation(float f) {
		  hpsPFTauDiscriminationByVLooseIsolation = f;
		}
		void SetHpsPFTauDiscriminationByLooseIsolation(float f) {
		  hpsPFTauDiscriminationByLooseIsolation = f;
		}
		void SetHpsPFTauDiscriminationByMediumIsolation(float f) {
		  hpsPFTauDiscriminationByMediumIsolation = f;
		}
		void SetHpsPFTauDiscriminationByTightIsolation(float f) {
		  hpsPFTauDiscriminationByTightIsolation = f;
		}

		void SetHpsPFTauDiscriminationByVLooseIsolationDBSumPtCorr(float f) {
		  hpsPFTauDiscriminationByVLooseCombinedIsolationDBSumPtCorr = f;
		}
		void SetHpsPFTauDiscriminationByLooseIsolationDBSumPtCorr(float f) {
		  hpsPFTauDiscriminationByLooseIsolationDBSumPtCorr = f;
		}
		void SetHpsPFTauDiscriminationByMediumIsolationDBSumPtCorr(float f) {
		  hpsPFTauDiscriminationByMediumIsolationDBSumPtCorr = f;
		}
		void SetHpsPFTauDiscriminationByTightIsolationDBSumPtCorr(float f) {
		  hpsPFTauDiscriminationByTightIsolationDBSumPtCorr = f;
		}


		void SetHpsPFTauDiscriminationByVLooseCombinedIsolationDBSumPtCorr(float f) {
		  hpsPFTauDiscriminationByVLooseCombinedIsolationDBSumPtCorr = f;
		}
		void SetHpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr(float f) {
		  hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr = f;
		}
		void SetHpsPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr(float f) {
		  hpsPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr = f;
		}
		void SetHpsPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr(float f) {
		  hpsPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr = f;
		}


		ClassDef(TCTau, 5);

};

#endif

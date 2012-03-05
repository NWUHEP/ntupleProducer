#ifndef _TGENPARTICLE_H
#define _TGENPARTICLE_H

#include "TObject.h"
#include "TLorentzVector.h"
#include "TVector2.h"
#include "TVector3.h"
#include <vector>

class TCGenParticle : public TLorentzVector {
	private:
		TVector3 _position;
		TLorentzVector _p4;
		int charge;
		int mother;
        int grandmother;
		int PDGID;
		//    std::vector<int> daughters;

	public:
		TCGenParticle();
		virtual ~TCGenParticle();

		TVector3 Position() const;
		TLorentzVector P4() const;
		TVector2 P2() const;
		float Et() const;
		float Pt() const;
		int Charge() const;
		int Mother();
        int Grandmother();
		int GetPDGId();
		//    std::vector<int> GetDaughters();

		// "set" methods ---------
		void SetPosition(float x, float y, float z);
		void SetP4(TLorentzVector p4);
		void SetP4(float px, float py, float pz, float e);
		void SetCharge(int c);
		//    void AddDaughter(int d);
		void SetMother(int m);
        void SetGrandmother(int g);
		void SetPDGId(int pdg_id);

		ClassDef(TCGenParticle, 1);

};

#endif	

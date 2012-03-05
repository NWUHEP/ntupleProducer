/* 
 * File:   TCGenJet.h
 * Author: Nate 0.
 *
 * Created on Nov 29, 2010, 1_39 PM
 */

#ifndef _TCGENJET_H
#define	_TCGENJET_H

#include "TObject.h"
#include "TLorentzVector.h"
#include "TVector2.h"

class TCGenJet : public TLorentzVector {
private:
    TLorentzVector _p4;
    TLorentzVector _progenitorP4;
    TVector3 _vtx;
    //    TVector3 _assocPV;

    float _hadEnergy;
    float _emEnergy;
    float _invEnergy;
    float _auxEnergy;

    unsigned int _numConstit;
    unsigned int _numChPart;
    int _jetFlavor;

public:
    TCGenJet();
    virtual ~TCGenJet();

    // "get" methods -----------

    TLorentzVector P4() const;
    TLorentzVector ProgenitorP4() const;
    TVector2 P2() const;
    float Et() const;
    float Pt() const;

    float HadEnergy() const;
    float EmEnergy() const;
    float InvEnergy() const;
    float AuxEnergy() const;    
    int   JetFlavor() const;
    unsigned int NumConstit() const;
    unsigned int NumChPart() const;

    TVector3 Vtx() const;
    //    TVector3 AssocVtx() const;

    // "set" methods ---------
    void SetP4(TLorentzVector p4);
    void SetP4(float px, float py, float pz, float e);
    void SetProgenitorP4(TLorentzVector p4);
    void SetVtx(float vx, float vy, float vz);
    //  void SetAssocVtx(float vx, float vy, float vz);

    void SetHadEnergy(float h);
    void SetEmEnergy(float e);
    void SetInvEnergy(float i);
    void SetAuxEnergy(float a);
    void SetNumConstit(unsigned int n);
    void SetNumChPart(unsigned int n);
    void SetJetFlavor(int f);

    ClassDef(TCGenJet, 1);

};

#endif	/* _TCGENJET_H */


/* 
 * File:   TCGenJet.h
 * Author: Nate 0oooo).
 *
 * Created on Nov 29, 2010, 1_39 PM
 */

#ifndef _TCGENJET_H
#define	_TCGENJET_H

#include <iostream>
#include "TObject.h"
#include "TLorentzVector.h"
#include "TCPhysObject.h"

class TCGenJet : public TCPhysObject {
private:
    TLorentzVector _progenitorP4;

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

    TLorentzVector ProgenitorP4() const;

    float HadEnergy() const;
    float EmEnergy() const;
    float InvEnergy() const;
    float AuxEnergy() const;    
    int   JetFlavor() const;
    unsigned int NumConstit() const;
    unsigned int NumChPart() const;

    // "set" methods ---------

    void SetProgenitorP4(TLorentzVector p4);

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


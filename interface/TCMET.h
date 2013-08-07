/* 
 * File:   TCMET.h
 * Author: Anton A.
 *
 * Created on April 30, 2010, 2:49 PM
 */

#ifndef _TCMET_H
#define	_TCMET_H

#include "TObject.h"
#include "TLorentzVector.h"
#include "TVector2.h"

class TCMET : public TVector2 {
private:

    TVector2 _genMET;

    float _sumEt;
    float _muonFraction;
    float _neutralHadronFraction;
    float _neutralEMFraction;
    float _chargedHadronFraction;
    float _chargedEMFraction;
    float _unCorPhi;

//Significance
    float _Significance;
    float _SigmaX2;

public:
    TCMET();
    virtual ~TCMET();

    // "get" methods -----------

    float SumEt() const;
    float MuonFraction() const;
    float NeutralHadronFraction() const;
    float NeutralEMFraction() const;
    float ChargedHadronFraction() const;
    float ChargedEMFraction() const;
    float UnCorPhi() const;

//Significance
    float Significance() const;
    float SigmaX2() const;

    // "set" methods ---------

    void SetSumEt(float n);
    void SetMuonFraction(float n);
    void SetNeutralHadronFraction(float n);
    void SetNeutralEMFraction(float n);
    void SetChargedHadronFraction(float n);
    void SetChargedEMFraction(float n);
    void SetUnCorPhi(float n);

//Significance
    void SetSignificance(float n);
    void SetSigmaX2(float n);

    ClassDef(TCMET, 1);

};

#endif	/* _TCMET_H */


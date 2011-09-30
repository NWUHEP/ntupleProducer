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

class TCMET : public TObject {
private:

    float _sumEt;
    float _met;
    float _phi;
    float _corSumEt;
    float _corMet;
    float _corPhi;
    float _photonEtFraction;
    float _electronEtFraction;
    float _muonEtFraction;
    float _neutralHadronEtFraction;
    float _chargedHadronEtFraction;
    float _hfEMEtFraction;
    float _hfHadronEtFraction;

public:
    TCMET();
    virtual ~TCMET();

    // "get" methods -----------

    float SumEt() const;
    float Met() const;
    float Phi() const;
    float CorrectedSumEt() const;
    float CorrectedMet() const;
    float CorrectedPhi() const;
    float PhotonEtFraction() const;
    float ElectronEtFraction() const;
    float MuonEtFraction() const;
    float NeutralHadronEtFraction() const;
    float ChargedHadronEtFraction() const;
    float HFEMEtFraction() const;
    float HFHadronEtFraction() const;
    float Significance() const;

    // "set" methods ---------

    void SetSumEt(float n);
    void SetMet(float n);
    void SetPhi(float n);
    void SetCorrectedSumEt(float n);
    void SetCorrectedMet(float n);
    void SetCorrectedPhi(float n);
    void SetPhotonEtFraction(float n);
    void SetElectronEtFraction(float n);
    void SetMuonEtFraction(float n);
    void SetNeutralHadronEtFraction(float n);
    void SetChargedHadronEtFraction(float n);
    void SetHFEMEtFraction(float n);
    void SetHFHadronEtFraction(float n);

    ClassDef(TCMET, 1);

};

#endif	/* _TCMET_H */


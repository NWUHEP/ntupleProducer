#ifndef _TCELECTRON_H
#define	_TCELECTRON_H

#include "TObject.h"
#include "TLorentzVector.h"

class TCElectron : public TObject {
private:
    TLorentzVector _p4;
    TVector3 _vtx;
    int _charge;
    float _hadOverEm;
    float _dPhiSuperCluster;
    float _dEtaSuperCluster;
    float _sigmaIetaIeta;

    float _emIso03;
    float _hadIso03;
    float _trkIso03;
    float _emIso04;
    float _hadIso04;
    float _trkIso04;

    float _ptError;

    float _convDist;
    float _convDcot;
    float _convRad;
    int _convFlag;

    bool _isEB;        // true if particle is in ECAL Barrel
    bool _isEE;        // true if particle is in ECAL Endcaps
    bool _isInGap;

    float _normalizedChi2;
    int _numberOfValidPixelHits;
    int _numberOfValidTrackerHits;
    int _numberOfLostPixelHits;
    int _numberOfLostTrackerHits;

    float _pfIso_Pt03;
    float _pfIso_Neutral03;
    float _pfIso_Gamma03;
    float _pfIso_Pt04;
    float _pfIso_Neutral04;
    float _pfIso_Gamma04;
    float _pfIso_Pt05;
    float _pfIso_Neutral05;
    float _pfIso_Gamma05;

    int _cut95;
    int _cut90;
    int _cut85;
    int _cut80;
    int _cut70;
    int _cut60;


public:
    TCElectron();
    virtual ~TCElectron();

    // "get" methods -----------

    TLorentzVector P4() const;
    TVector2 P2() const;
    float Pt() const;
    float PtError() const;
    TVector3 Vtx() const;
    float Et() const;
    float Eta() const;
    float Phi() const;
    int Charge() const;
    float EmIso() const;
    float HadIso() const;
    float TrkIso() const;
    float EmIso03() const;
    float HadIso03() const;
    float TrkIso03() const;
    float EmIso04() const;
    float HadIso04() const;
    float TrkIso04() const;
    float HadOverEm() const;
    float DphiSuperCluster() const;
    float DetaSuperCluster() const;
    float SigmaIetaIeta() const;

    float NormalizedChi2() const;

    float PfRelIso(float coneSize) const;
    float PfSumPt(float coneSize) const;
    float PfEGamma(float coneSize) const;
    float PfENeutral(float coneSize) const;


    int ConversionFlag() const;
    float ConversionDist() const;
    float ConversionDcot() const;
    float ConversionRad() const;

    bool IsEB() const;
    bool IsEE() const;
    bool IsInGap() const;

    int NumberOfValidPixelHits() const;
    int NumberOfValidTrackerHits() const;
    int NumberOfLostPixelHits() const;
    int NumberOfLostTrackerHits() const;

    int CutLevel(int lvl) const;
    float Dxy(TVector3 *primVtx) const;
    float Dz(TVector3 *primVtx) const;

    //--------------------------
    // "set" methods ---------
    //--------------------------

    void SetP4(TLorentzVector p4);
    void SetP4(float px, float py, float pz, float e);
    void SetVtx(float vx, float vy, float vz);

    void SetCharge(int c);

    void SetEmIso03(float e);
    void SetHadIso03(float h);
    void SetTrkIso03(float t);
    void SetEmIso04(float e);
    void SetHadIso04(float h);
    void SetTrkIso04(float t);
 
    void SetHadOverEm(float h);
    void SetDphiSuperCluster(float dp);
    void SetDetaSuperCluster(float de);
    void SetSigmaIetaIeta(float sieie);
    void SetConversionDist(float d);
    void SetConversionDcot(float d);
    void SetConversionRad(float r);
    void SetConversionFlag(int f);

    void SetNumberOfValidPixelHits(int n);
    void SetNumberOfValidTrackerHits(int n);
    void SetNumberOfLostPixelHits(int n);
    void SetNumberOfLostTrackerHits(int n);
    void SetNormalizedChi2(float n);

    void SetPfSumPt(float coneSize, float f); 
    void SetPfEGamma(float coneSize, float f);
    void SetPfENeutral(float coneSize, float f);

    void SetIsEB(bool b);
    void SetIsEE(bool b);
    void SetIsInGap(bool b);

    void SetCutLevel(int cut, int lvl);

ClassDef(TCElectron, 1);

};

#endif	/* _TCELECTRON_H */



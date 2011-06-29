#ifndef _TCELECTRON_H
#define	_TCELECTRON_H

#include "TObject.h"
#include "TLorentzVector.h"

class TCElectron : public TObject {
private:
    TLorentzVector _p4;
    TVector3 _vtx;
    int _charge;
    float _emIso;
    float _hadIso;
    float _trkIso;
    float _hadOverEm;
    float _dPhiSuperCluster;
    float _dEtaSuperCluster;
    float _sigmaIetaIeta;

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
    float pt() const;
    float ptError() const;
    TVector3 Vtx() const;
    float Et() const;
    float eta() const;
    float phi() const;
    int charge() const;
    float emIso() const;
    float hadIso() const;
    float trkIso() const;
    float hadOverEm() const;
    float dPhiSuperCluster() const;
    float dEtaSuperCluster() const;
    float sigmaIetaIeta() const;

    float normalizedChi2() const;

    float pfRelIso(float coneSize) const;
    float pfSumPt(float coneSize) const;
    float pfEGamma(float coneSize) const;
    float pfENeutral(float coneSize) const;


    int conversionFlag() const;
    float conversionDist() const;
    float conversionDcot() const;
    float conversionRad() const;

    bool isEB() const;
    bool isEE() const;
    bool isInGap() const;

    int numberOfValidPixelHits() const;
    int numberOfValidTrackerHits() const;
    int numberOfLostPixelHits() const;
    int numberOfLostTrackerHits() const;

    int CutLevel(int lvl) const;
    // float dxy(TVector3 *primVtx) const;
    // float dz(TVector3 *primVtx) const;

    //--------------------------
    // "set" methods ---------
    //--------------------------

    void SetP4(TLorentzVector p4);
    void SetP4(float px, float py, float pz, float e);
    void SetVtx(float vx, float vy, float vz);

    void SetCharge(int c);
    void SetEMIso(float e);
    void SetHADIso(float h);
    void SetTRKIso(float t);
 
    void SetHadOverEm(float h);
    void SetDPhiSuperCluster(float dp);
    void SetDEtaSuperCluster(float de);
    void SetSigmaIetaIeta(float sieie);
    void SetConversionDist(float d);
    void SetConversionDcot(float d);
    void SetConversionRad(float r);
    void SetConversionFlag(int f);

    void SetnumberOfValidPixelHits(int n);
    void SetnumberOfValidTrackerHits(int n);
    void SetnumberOfLostPixelHits(int n);
    void SetnumberOfLostTrackerHits(int n);
    void SetnormalizedChi2(float n);

    void SetPFSumPt(float coneSize, float f); 
    void SetPFEGamma(float coneSize, float f);
    void SetPFENeutral(float coneSize, float f);

    void SetisEB(bool b);
    void SetisEE(bool b);
    void SetisInGap(bool b);

    void SetCutLevel(int cut, int lvl);

ClassDef(TCElectron, 1);

};

#endif	/* _TCELECTRON_H */



#ifndef _TCMUON_H
#define	_TCMUON_H

#include "TObject.h"
#include "TLorentzVector.h"
#include "TVector2.h"
#include "TVector3.h"

class TCMuon : public TObject {
private:
    TLorentzVector _p4;
    TVector3 _vtx;
    int _charge;
    bool _isTRK;
    bool _isGLB;
    float _caloComp;
    float _segComp;
    float _emIso;
    float _hadIso;
    float _trkIso;
    float _ptError; 

    int _numberOfMatches;
    int _nTracks;
    int _numberOfValidPixelHits;
    int _numberOfValidTrackerHits;
    int _numberOfLostPixelHits;
    int _numberOfLostTrackerHits;
    int _numberOfValidMuonHits;
    float _normalizedChi2;

    float _pfIso_Pt03;
    float _pfIso_Neutral03;
    float _pfIso_Gamma03;
    float _pfIso_Pt04;
    float _pfIso_Neutral04;
    float _pfIso_Gamma04;
    float _pfIso_Pt05;
    float _pfIso_Neutral05;
    float _pfIso_Gamma05;

public:
    TCMuon();
    virtual ~TCMuon();

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
    bool isGLB() const;
    bool isTRK() const;
    float caloComp() const;
    float segComp() const;
    float emIso() const;
    float hadIso() const;
    float trkIso() const;

    int nTracks() const;
    int numberOfValidPixelHits() const;
    int numberOfValidTrackerHits() const;
    int numberOfLostPixelHits() const;
    int numberOfLostTrackerHits() const;
    int numberOfValidMuonHits() const;
    float normalizedChi2() const;
    int numberOfMatches() const;

    float pfRelIso(float coneSize) const;
    float pfSumPt(float coneSize) const;
    float pfEGamma(float coneSize) const;
    float pfENeutral(float coneSize) const;

  // float dxy(TVector3 *primVtx) const;
  // float dz(TVector3 *primVtx) const;

   // "set" methods ---------
    void SetP4(TLorentzVector p4);
    void SetP4(float px, float py, float pz, float e);
    void SetptError(float er);
    void SetVtx(float vx, float vy, float vz);

    void SetnTracks(int n);
    void SetCharge(int c);
    void SetisGLB(bool t);
    void SetisTRK(bool t);
    void SetCaloComp(float c);
    void SetSegComp(float s);
    void SetEMIso(float e);
    void SetHADIso(float h);
    void SetTRKIso(float t);
    void SetnumberOfMatches(int n);
    void SetnumberOfValidPixelHits(int n);
    void SetnumberOfValidTrackerHits(int n);
    void SetnumberOfValidMuonHits(int n);
    void SetnumberOfLostPixelHits(int n);
    void SetnumberOfLostTrackerHits(int n);
    void SetnormalizedChi2(float n);

    void SetPFSumPt(float coneSize, float f); 
    void SetPFEGamma(float coneSize, float f);
    void SetPFENeutral(float coneSize, float f);

    ClassDef(TCMuon, 1);

};

#endif	/* _TCMUON_H */



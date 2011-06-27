#ifndef _NUMUON_H
#define	_NUMUON_H

#include "TObject.h"
#include "TLorentzVector.h"

class TCMuon : public TObject {
private:
    TLorentzVector _p4;
    TVector3 _vtx;
    //    TVector3 _assocPV;
    float _eta;
    float _phi;
    int _charge;
    bool _isTRK;
    bool _isGLB;
    float _dxy;
    int _nPXLHits;
    int _nTRKHits;
    int _nMatchSeg;
    int _nValidMuHits;
    float _normChi2;
    float _caloComp;
    float _segComp;
    float _emIso;
    float _hadIso;
    float _trkIso;

public:
    TCMuon();
    virtual ~TCMuon();

    // "get" methods -----------

    TLorentzVector p4() const;
    float pT() const;
    TVector3 Vtx() const;
    float eta() const;
    float phi() const;
    int charge() const;
    bool isGLB() const;
    bool isTRK() const;
    float dxy() const;
    int nPXLHits() const;
    int nTRKHits() const;
    int nMatchSeg() const;
    int nValidMuHits() const;
    float normChi2() const;
    float caloComp() const;
    float segComp() const;
    float emIso() const;
    float hadIso() const;
    float trkIso() const;
    //    TVector3 AssocVtx() const;

    // "set" methods ---------
    void Setp4(TLorentzVector p4);
    void Setp4(float px, float py, float pz, float e);
    void SetVtx(float vx, float vy, float vz);
    //  void SetAssocVtx(float vx, float vy, float vz);

    void SetEta(float e);
    void SetPhi(float p);
    void SetCharge(int c);
    void SetisGLB(bool t);
    void SetisTRK(bool t);
    void Setdxy(float d);
    void SetnPXLHits(int n);
    void SetnTRKHits(int n);
    void SetnValidMuHits(int n);
    void SetnMatchSeg(int n);
    void SetNormChi2(float c);
    void SetCaloComp(float c);
    void SetSegComp(float s);
    void SetEMIso(float e);
    void SetHADIso(float h);
    void SetTRKIso(float t);

    ClassDef(TCMuon, 1);

};

#endif	/* _NUMUON_H */



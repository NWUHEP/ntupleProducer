#ifndef _TCPHOTON_H
#define	_TCPHOTON_H

#include "TObject.h"
#include "TLorentzVector.h"

class TCPhoton : public TObject {
private:
    TLorentzVector _p4;
    TVector3 _vtx;
    int   _charge;
    float _normChi2;
    float _emIso; // 
    float _hadIso; // 
    float _trkIso; // 
    float _hadOverEm; // 
    float _dPhiSuperCluster;
    float _dEtaSuperCluster;
    float _sigmaIEtaIEta; // 
    float _r9;
    float _sigmaIPhiIPhi; 
	float _e2OverE9;
	float _etaSupercluster;
    bool  _trackVeto;

public:
    TCPhoton();
    virtual ~TCPhoton();

    // "get" methods -----------

    TLorentzVector P4() const;
    TVector3 Vtx() const;
    float Pt() const;
    float Eta() const;
    float Phi() const;
    int   Charge() const;
    float NormChi2() const;
    float EmIso() const;
    float HadIso() const;
    float TrkIso() const;
    float HadOverEm() const;
    float DPhiSuperCluster() const;
    float DEtaSuperCluster() const;
    float SigmaIEtaIEta() const;
    float SigmaIPhiIPhi() const;
    float R9() const; 
    float E2OverE9() const; 
    float EtaSupercluster() const;
    bool  TrackVeto() const;

    //    TVector3 AssocVtx() const;

    // "set" methods ---------
    void SetP4(TLorentzVector p4);
    void SetP4(float px, float py, float pz, float e);
    void SetVtx(float vx, float vy, float vz);
    //  void SetAssocVtx(float vx, float vy, float vz);

    void SetCharge(int c);
    void SetNormChi2(float c);
    void SetEMIso(float e);
    void SetHADIso(float h);
    void SetTRKIso(float t);
 
    void SetHadOverEm(float h);
    void SetDPhiSuperCluster(float d);
    void SetDEtaSuperCluster(float d);
    void SetSigmaIEtaIEta(float s);
    void SetR9(float r);
    void SetSigmaIPhiIPhi(float s);
    void Sete2OverE9(float e);
    void SetEtaSupercluster(float e);
    void SetTrackVeto(bool t);

    ClassDef(TCPhoton, 1);

};

#endif



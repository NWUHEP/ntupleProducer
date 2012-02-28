#ifndef _TCPHOTON_H
#define	_TCPHOTON_H

#include "TObject.h"
#include "TLorentzVector.h"
#include "TArrayF.h"

class TCPhoton : public TLorentzVector {
private:
    TLorentzVector _p4;
    //std::pair<TLorentzVector, TLorentzVector> _convP4;
    TVector3 _vtx;
    TArrayF  _vtxIso;

    int   _charge;
    float _normChi2;
    float _emIso; // 
    float _hadIso; // 
    float _trkIso; // 
    float _pfIsoNeutral;
    float _pfIsoCharged;
    float _pfIsoPhoton;
    float _hadOverEm; // 
    float _dPhiSuperCluster;
    float _dEtaSuperCluster;
    float _sigmaIEtaIEta; // 
    float _r9;
    float _sigmaIPhiIPhi; 
	float _e2OverE9;
	float _etaSupercluster;
    bool  _trackVeto;

    //conversion info
    int   _nConversions; //
    float _conversionDz; //
    float _conversionDxy; //



public:
    TCPhoton();
    TCPhoton(int size);
    virtual ~TCPhoton();

    // "get" methods -----------

    TLorentzVector P4() const;
    TVector3 Vtx() const;
    float Pt() const;
    float Eta() const;
    float Phi() const;

    int   Charge() const;
    float NormChi2() const;
    bool  TrackVeto() const;

    float EmIso() const;
    float HadIso() const;
    float TrkIso() const;
    float PFIsoNeutral() const;
    float PFIsoCharged() const;
    float PFIsoPhoton() const;
    float VertexIso(int iVtx) const;

    float HadOverEm() const;
    float DPhiSuperCluster() const;
    float DEtaSuperCluster() const;
    float SigmaIEtaIEta() const;
    float SigmaIPhiIPhi() const;
    float R9() const; 
    float E2OverE9() const; 
    float EtaSupercluster() const;

    int    NumberOfConversions() const;
    float  ConversionDz() const;
    float  ConversionDxy() const;
    //std::pair<TLorentzVector, TLorentzVector>  ConversionPairP4() const;

    //    TVector3 AssocVtx() const;

    // "set" methods ---------
    void SetP4(TLorentzVector p4);
    void SetP4(float px, float py, float pz, float e);
    void SetVtx(float vx, float vy, float vz);

    void SetCharge(int c);
    void SetNormChi2(float c);
    void SetEMIso(float e);
    void SetHADIso(float h);
    void SetTRKIso(float t);
    void SetPFIsoNeutral(float n);
    void SetPFIsoCharged(float c);
    void SetPFIsoPhoton(float p);
 
    void SetHadOverEm(float h);
    void SetDPhiSuperCluster(float d);
    void SetDEtaSuperCluster(float d);
    void SetSigmaIEtaIEta(float s);
    void SetR9(float r);
    void SetSigmaIPhiIPhi(float s);
    void Sete2OverE9(float e);
    void SetEtaSupercluster(float e);
    void SetTrackVeto(bool t);

    void SetNumberOfConversions(int n);
    void SetConversionDz(float d);
    void SetConversionDxy(float d);
    //void SetConversionPairP4(TLorentzVector p1, TLorentzVector p2);

    void SetVertexIso(int iVtx, float iso);

    ClassDef(TCPhoton, 1);
};

#endif



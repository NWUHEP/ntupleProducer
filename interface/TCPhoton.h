#ifndef _TCPHOTON_H
#define	_TCPHOTON_H

#include "TObject.h"
#include "TLorentzVector.h"
#include "TArrayF.h"
#include <vector>

using namespace std;

class TCPhoton : public TObject {
private:
    TLorentzVector _p4;
    //std::pair<TLorentzVector, TLorentzVector> _convP4;
    TVector3 _vtx;
    int   _charge;

    float _customEm;
    float _customNh;
    float _customCh;

    float _normChi2;
    float _emIsoDR04; // 
    float _hadIsoDR04; // 
    float _trkIsoDR04; // 
    vector<float> _trkIsoVtxDR04;
    float _emIsoDR03; // 
    float _hadIsoDR03; // 
    float _trkIsoDR03; // 
    vector<float> _trkIsoVtxDR03;
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
    virtual ~TCPhoton();

    // "get" methods -----------

    float CustomEm() const;
    float CustomNh() const;
    float CustomCh() const;

    TLorentzVector P4() const;
    TVector3 Vtx() const;
    float Pt() const;
    float Eta() const;
    float Phi() const;
    int   Charge() const;
    float NormChi2() const;
    float EmIsoDR04() const;
    float HadIsoDR04() const;
    float TrkIsoDR04() const;
    vector<float> TrkIsoVtxDR04() const;
    float EmIsoDR03() const;
    float HadIsoDR03() const;
    float TrkIsoDR03() const;
    vector<float> TrkIsoVtxDR03() const;
    float PFIsoNeutral() const;
    float PFIsoCharged() const;
    float PFIsoPhoton() const;
    float HadOverEm() const;
    float DPhiSuperCluster() const;
    float DEtaSuperCluster() const;
    float SigmaIEtaIEta() const;
    float SigmaIPhiIPhi() const;
    float R9() const; 
    float E2OverE9() const; 
    float EtaSupercluster() const;
    bool  TrackVeto() const;

    int    NumberOfConversions() const;
    float  ConversionDz() const;
    float  ConversionDxy() const;
    //std::pair<TLorentzVector, TLorentzVector>  ConversionPairP4() const;

    //    TVector3 AssocVtx() const;

    // "set" methods ---------
    void SetP4(TLorentzVector p4);
    void SetP4(float px, float py, float pz, float e);
    void SetVtx(float vx, float vy, float vz);

    void SetCustomIso(float em, float nh, float ch);

    void SetCharge(int c);
    void SetNormChi2(float c);
    void SetEMIsoDR04(float e);
    void SetHADIsoDR04(float h);
    void SetTRKIsoDR04(float t);
    void SetTRKIsoVtxDR04(float t);
    void SetEMIsoDR03(float e);
    void SetHADIsoDR03(float h);
    void SetTRKIsoDR03(float t);
    void SetTRKIsoVtxDR03(float t);
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

    ClassDef(TCPhoton, 1);
};

#endif



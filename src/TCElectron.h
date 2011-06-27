#ifndef _NUELECTRON_H
#define	_NUELECTRON_H

#include "TObject.h"
#include "TLorentzVector.h"

class TCElectron : public TObject {
private:
    TLorentzVector _p4;
    TVector3 _vtx;
    //    TVector3 _assocPV;
    int _charge;
	 int _idMap;
	 int _convFlag;
    float _dxy;
    float _normChi2;
    float _emIso;
    float _hadIso;
    float _trkIso;
    float _HoverE;
    float _dPhiSC;
	 float _dEtaSC;
	 float _sig_IEtaIEta;
	 float _pfChargedHadronIso;
	 float _pfNeutralHadronIso;
	 float _pfPhotonIso;
	 float _convDist;
	 float _convDcot;
	 float _convRad;

public:
    TCElectron();
    virtual ~TCElectron();

    // "get" methods -----------

    TLorentzVector p4() const;
    float pT() const;
    TVector3 Vtx() const;
    int charge() const;
	 int idMap() const;
	 int conversionFlag() const;
    float dxy() const;
    float eta() const;
    float phi() const;
    float normChi2() const;
    float emIso() const;
    float hadIso() const;
    float trkIso() const;
    float HoverE() const;
    float dPhiSC() const;
    float dEtaSC() const;
    float Sig_IEtaIEta() const;
    float pfChargedHadronIso() const;
    float pfNeutralHadronIso() const;
    float pfPhotonIso() const;
	 float conversionDist() const;
	 float conversionDcot() const;
	 float conversionRad() const;

    //    TVector3 AssocVtx() const;

    // "set" methods ---------
    void Setp4(TLorentzVector p4);
    void Setp4(float px, float py, float pz, float e);
    void SetVtx(float vx, float vy, float vz);
    //  void SetAssocVtx(float vx, float vy, float vz);

    void SetCharge(int c);
	 void SetIDMap(int i);
	 void SetConversionFlag(int f);
    void Setdxy(float d);
    void SetNormChi2(float c);
    void SetEMIso(float e);
    void SetHADIso(float h);
    void SetTRKIso(float t);
    void SetHoverE(float h);
    void SetdPhiSC(float d);
    void SetdEtaSC(float d);
    void SetSig_IEtaIEta(float s);
    void SetPFChargedHadronIso(float c);
    void SetPFNeutralHadronIso(float n);
    void SetPFPhotonIso(float g);
	 void SetConversionDist(float d);
	 void SetConversionDcot(float d);
	 void SetConversionRad(float r);

    ClassDef(TCElectron, 1);

};

#endif	/* _NUELECTRON_H */



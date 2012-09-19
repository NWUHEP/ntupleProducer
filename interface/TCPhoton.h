#ifndef _TCPHOTON_H
#define	_TCPHOTON_H

#include "TObject.h"
#include "TLorentzVector.h"
#include "TArrayF.h"
#include "TCPhysObject.h"
#include <vector>

using namespace std;

class TCPhoton : public TCPhysObject {
private:

    // ID variables
    float _normChi2;
    float _hadOverEm; // 
    float _sigmaIEtaIEta; // 
    float _r9;
    float _sigmaIPhiIPhi; 
    float _e2OverE9;
    bool  _trackVeto;

    // supercluster information
    float _dPhiSC;
    float _dEtaSC;
    float _etaSC;
    float _energySC;

    //conversion info
    int   _nConversions; //
    float _conversionDz; //
    float _conversionDxy; //

    // vertex-by-vertex iso
    vector<float> _trkIsoVtxDR03;
    vector<float> _trkIsoVtxDR04;

public:
    TCPhoton();
    virtual ~TCPhoton();

    // "get" methods -----------

    float NormChi2() const;
    float HadOverEm() const;
    float SigmaIEtaIEta() const;
    float SigmaIPhiIPhi() const;
    float R9() const; 
    float E2OverE9() const; 
    bool  TrackVeto() const;

    float DPhiSC() const;
    float DEtaSC() const;
    float EnergySC() const;
    float EtaSC() const;

    int    NumberOfConversions() const;
    float  ConversionDz() const;
    float  ConversionDxy() const;

    float TrkIsoVtxDR03(int);
    float TrkIsoVtxDR04(int);

    // "set" methods ---------

    void SetNormChi2(float c);
    void SetHadOverEm(float h);
    void SetSigmaIEtaIEta(float s);
    void SetSigmaIPhiIPhi(float s);
    void SetR9(float r);
    void SetE2OverE9(float e);
    void SetTrackVeto(bool t);

    void SetDPhiSC(float d);
    void SetDEtaSC(float d);
    void SetEtaSC(float e);
    void SetEnergySC(float e);

    void SetNumberOfConversions(int n);
    void SetConversionDz(float d);
    void SetConversionDxy(float d);

    void SetTRKIsoVtxDR04(float);
    void SetTRKIsoVtxDR03(float);
 
    ClassDef(TCPhoton, 1);
};

#endif



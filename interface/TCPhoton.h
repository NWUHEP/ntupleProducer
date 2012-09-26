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
    bool    _convVeto;

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

    bool  ConversionVeto() const;

    // "set" methods ---------

    void SetNormChi2(float);
    void SetHadOverEm(float);
    void SetSigmaIEtaIEta(float);
    void SetSigmaIPhiIPhi(float);
    void SetR9(float);
    void SetE2OverE9(float);
    void SetTrackVeto(bool);

    void SetDPhiSC(float);
    void SetDEtaSC(float);
    void SetEtaSC(float);
    void SetEnergySC(float);

    void SetConversionVeto(bool);

    ClassDef(TCPhoton, 1);
};

#endif



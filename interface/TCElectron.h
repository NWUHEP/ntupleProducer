#ifndef _TCELECTRON_H
#define	_TCELECTRON_H

#include "TObject.h"
#include "TLorentzVector.h"
#include "TCPhysObject.h"

class TCElectron : public TCPhysObject {
    private:

        float _ptError;
        float _hadOverEm;
        float _dPhiSuperCluster;
        float _dEtaSuperCluster;
        float _sigmaIetaIeta;
        float _eOverP;
        float _fBrem;
        float _r9;

        float _scEta;

        bool  _convVeto;
        short _convMissHits;

        bool _isEB;        // true if particle is in ECAL Barrel
        bool _isEE;        // true if particle is in ECAL Endcaps
        bool _isInGap;

        float _normalizedChi2;
        int _numberOfValidPixelHits;
        int _numberOfValidTrackerHits;
        int _numberOfLostPixelHits;
        int _numberOfLostTrackerHits;

        int _cut95;
        int _cut90;
        int _cut85;
        int _cut80;
        int _cut70;
        int _cut60;

        TLorentzVector _regressionMomCombP4;


    public:
        TCElectron();
        virtual ~TCElectron();

        // "get" methods -----------

        float PtError() const;

        float HadOverEm() const;
        float DphiSuperCluster() const;
        float DetaSuperCluster() const;
        float SigmaIEtaIEta() const;
        float FBrem() const;
        float EOverP() const;
        float NormalizedChi2() const;

        float SCEta() const;
        float R9() const; 

        bool  ConversionVeto() const;
        short ConversionMissHits() const;

        bool IsEB() const;
        bool IsEE() const;
        bool IsInGap() const;

        int NumberOfValidPixelHits() const;
        int NumberOfValidTrackerHits() const;
        int NumberOfLostPixelHits() const;
        int NumberOfLostTrackerHits() const;

        int CutLevel(int lvl) const;
        bool PassID(int lvl) const;
        bool PassConversion(int lvl) const;
        bool PassIsolation(int lvl) const;

        TLorentzVector RegressionMomCombP4() const;

        //--------------------------
        // "set" methods ---------
        //--------------------------

        void SetPtError(float e);
        void SetHadOverEm(float h);
        void SetDphiSuperCluster(float dp);
        void SetDetaSuperCluster(float de);
        void SetSigmaIEtaIEta(float sieie);
        void SetEOverP(float e);
        void SetFBrem(float fb);

        void SetSCEta(float);

        void SetConversionVeto(bool);
        void SetConversionMissHits(short);

        void SetNumberOfValidPixelHits(int n);
        void SetNumberOfValidTrackerHits(int n);
        void SetNumberOfLostPixelHits(int n);
        void SetNumberOfLostTrackerHits(int n);
        void SetNormalizedChi2(float n);
        void SetR9(float r);

        void SetIsEB(bool b);
        void SetIsEE(bool b);
        void SetIsInGap(bool b);

        void SetCutLevel(int cut, int lvl);

        void SetRegressionMomCombP4(TLorentzVector tmpP4);

        ClassDef(TCElectron, 1);
};

#endif	/* _TCELECTRON_H */



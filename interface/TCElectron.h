#ifndef _TCELECTRON_H
#define	_TCELECTRON_H

#include "TObject.h"
#include "TLorentzVector.h"
#include "Higgs/ntupleProducer/interface/TCPhysObject.h"

class TCElectron : public TCPhysObject {
    private:

        float _ptError;
        float _hadOverEm;
        float _dPhiSuperCluster;
        float _dEtaSuperCluster;
        float _sigmaIetaIeta;
        float _eOverP;
        float _fBrem;


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

        float PtError() const;

        float HadOverEm() const;
        float DphiSuperCluster() const;
        float DetaSuperCluster() const;
        float SigmaIetaIeta() const;
        float FBrem() const;
        float EOverP() const;
        float NormalizedChi2() const;

        int ConversionFlag() const;
        float ConversionDist() const;
        float ConversionDcot() const;
        float ConversionRad() const;

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

        //--------------------------
        // "set" methods ---------
        //--------------------------

        void SetHadOverEm(float h);
        void SetDphiSuperCluster(float dp);
        void SetDetaSuperCluster(float de);
        void SetSigmaIetaIeta(float sieie);
        void SetEOverP(float e);
        void SetFBrem(float fb);

        void SetConversionDist(float d);
        void SetConversionDcot(float d);
        void SetConversionRad(float r);
        void SetConversionFlag(int f);

        void SetNumberOfValidPixelHits(int n);
        void SetNumberOfValidTrackerHits(int n);
        void SetNumberOfLostPixelHits(int n);
        void SetNumberOfLostTrackerHits(int n);
        void SetNormalizedChi2(float n);

        void SetIsEB(bool b);
        void SetIsEE(bool b);
        void SetIsInGap(bool b);

        void SetCutLevel(int cut, int lvl);

        ClassDef(TCElectron, 1);

};

#endif	/* _TCELECTRON_H */



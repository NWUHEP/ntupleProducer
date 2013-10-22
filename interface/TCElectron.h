#ifndef _TCELECTRON_H
#define	_TCELECTRON_H

#include "TObject.h"
#include "TLorentzVector.h"
#include "TCEGamma.h"

class TCElectron : public TCEGamma {
    private:

        float _ptError;

        float _fBrem;

        float _mvaID;
        float _regEne;
        float _regErr;

        bool  _convVeto;
        short _convMissHits;

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


        float FBrem() const;

        float NormalizedChi2() const;

        float MvaID() const; 
        float EnergyRegression() const; 
        float EnergyRegressionErr() const; 

        bool  ConversionVeto() const;
        short ConversionMissHits() const;


        int NumberOfValidPixelHits() const;
        int NumberOfValidTrackerHits() const;
        int NumberOfLostPixelHits() const;
        int NumberOfLostTrackerHits() const;

        int  CutLevel(int lvl) const;
        bool PassID(int lvl) const;
        bool PassConversion(int lvl) const;
        bool PassIsolation(int lvl) const;

        TLorentzVector RegressionMomCombP4() const;

        //--------------------------
        // "set" methods ---------
        //--------------------------

        void SetPtError(float e);
        
        void SetFBrem(float fb);

        void SetConversionVeto(bool);
        void SetConversionMissHits(short);

        void SetNumberOfValidPixelHits(int n);
        void SetNumberOfValidTrackerHits(int n);
        void SetNumberOfLostPixelHits(int n);
        void SetNumberOfLostTrackerHits(int n);
        void SetNormalizedChi2(float n);
        
        void SetMvaID(float m);
        void SetEnergyRegression(float e);
        void SetEnergyRegressionErr(float e);

        void SetCutLevel(int cut, int lvl);

        void SetRegressionMomCombP4(TLorentzVector tmpP4);

        ClassDef(TCElectron, 1);
};

#endif	/* _TCELECTRON_H */



#ifndef _TCMUON_H
#define	_TCMUON_H

#include "TObject.h"
#include "TLorentzVector.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TCPhysObject.h"

class TCMuon : public TCPhysObject {
    private:

        float _ptError; 
        bool _isTRK;
        bool _isGLB;
        bool _isSoft;
        bool _isTight;
        bool _isGood;
        bool _isGoodLoose;
        float _caloComp;
        float _segComp;

        int _numberOfMatches;
        int _numberOfMatchedStations;
        int _numberOfValidPixelHits;
        int _numberOfValidTrackerHits;
        int _numberOfLostPixelHits;
        int _numberOfLostTrackerHits;
        int _numberOfValidMuonHits;
        int _trackLayersWithMeasurement;
        int _pixelLayersWithMeasurement;
        float _normalizedChi2;
        float _normalizedChi2_tracker;


        float _pfIsoPU;
        float _pfIsoChargedHad;
        float _pfIsoChargedPart;
        float _pfIsoNeutral;
        float _pfIsoPhoton;


    public:
        TCMuon();
        virtual ~TCMuon();

        // "get" methods -----------

        float PtError() const;

        bool IsGLB() const;
        bool IsTRK() const;
        bool IsSoft() const;
        bool IsTight() const;
        bool IsGood() const;
        bool IsGoodLoose() const;
        float CaloComp() const;
        float SegComp() const;

        int NumberOfValidPixelHits() const;
        int NumberOfValidTrackerHits() const;
        int NumberOfLostPixelHits() const;
        int NumberOfLostTrackerHits() const;
        int NumberOfValidMuonHits() const;
        float NormalizedChi2() const;
        float NormalizedChi2_tracker() const;
        int NumberOfMatches() const;
        int NumberOfMatchedStations() const;
        int TrackLayersWithMeasurement() const;
        int PixelLayersWithMeasurement() const;

        float PfIsoPU() const;
        float PfIsoChargedPart() const;
        float PfIsoChargedHad() const;
        float PfIsoCharged() const; //this will return Had
        float PfIsoNeutral() const;
        float PfIsoPhoton() const;

        void SetPtError(float er);
        void SetIsGLB(bool t);
        void SetIsTRK(bool t);
        void SetIsSoft(bool t);
        void SetIsTight(bool t);
        void SetIsGood(bool g);
        void SetIsGoodLoose(bool g);
        void SetCaloComp(float c);
        void SetSegComp(float s);

        void SetNumberOfMatches(int n);
        void SetNumberOfMatchedStations(int n);
        void SetNumberOfValidPixelHits(int n);
        void SetNumberOfValidTrackerHits(int n);
        void SetNumberOfValidMuonHits(int n);
        void SetNumberOfLostPixelHits(int n);
        void SetNumberOfLostTrackerHits(int n);
        void SetNormalizedChi2(float n);
        void SetNormalizedChi2_tracker(float n);
        void SetTrackLayersWithMeasurement(int n);
        void SetPixelLayersWithMeasurement(int n);


        void SetPfIsoPU(float f);
        void SetPfIsoChargedHad(float f);
        void SetPfIsoChargedPart(float f);
        void SetPfIsoNeutral(float f);
        void SetPfIsoPhoton(float f);

        ClassDef(TCMuon, 1);
};

#endif	/* _TCMUON_H */



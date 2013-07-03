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
        bool _isPF;
        bool _isTRK;
        bool _isGLB;
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
        float _normalizedChi2;

    public:
        TCMuon();
        virtual ~TCMuon();

        // "get" methods -----------

        float PtError() const;

        bool IsGLB() const;
        bool IsTRK() const;
        bool IsPF() const;
        float CaloComp() const;
        float SegComp() const;

        int NumberOfValidPixelHits() const;
        int NumberOfValidTrackerHits() const;
        int NumberOfLostPixelHits() const;
        int NumberOfLostTrackerHits() const;
        int NumberOfValidMuonHits() const;
        float NormalizedChi2() const;
        int NumberOfMatches() const;
        int NumberOfMatchedStations() const;
        int TrackLayersWithMeasurement() const;

        void SetPtError(float er);
        void SetIsGLB(bool t);
        void SetIsTRK(bool t);
        void SetIsPF(bool t);
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
        void SetTrackLayersWithMeasurement(int n);

        ClassDef(TCMuon, 1);
};

#endif	/* _TCMUON_H */



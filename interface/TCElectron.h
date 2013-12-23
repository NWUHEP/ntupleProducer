#ifndef _TCELECTRON_H
#define	_TCELECTRON_H

#include "TObject.h"
#include "TLorentzVector.h"
#include "TCEGamma.h"
#include "TCTrack.h"

class TCElectron : public TCEGamma {
 public:
  class Track : public TCTrack{
                       ClassDef(Track, 1);
                     };

 private:

  float _normChi2kf;
  float _ptError;

  float _fBrem;
  float _inverseEnergyMomentumDiff;

  float _EoP;
  float _EoPout;
  float _ESeedOverP;

  float _ip3d;
  float _ip3dSig;

  float _deltaEtaSeedCluster;
  float _deltaPhiSeedCluster;

  float _mvaID_Old;
  float _mvaID_HZZ;
  float _regEne;
  float _regErr;

  bool  _passConvVeto;
  short _convMissHits;

  float _convDcot;
  float _convDist;
  float _convRadius;

  int _trackerLayersWithMeasurement;
  int _numberOfValidHits;
  int _numberOfValidPixelHits;
  int _numberOfValidTrackerHits;
  int _numberOfLostPixelHits;
  int _numberOfLostTrackerHits;


  TLorentzVector _regressionMomCombP4;
  float _effArea;

  vector<TCElectron::Track> _tracks;

 public:
  TCElectron();
  virtual ~TCElectron();

  // "get" methods -----------
  vector<TCElectron::Track> GetTracks() const;

  float NormalizedChi2() const;
  float NormalizedChi2Gsf() const;
  float NormalizedChi2Kf() const;

  float InverseEnergyMomentumDiff() const;

  float EoP() const;
  float EoPout() const;
  float ESeedOverP() const;

  float IP3d() const;
  float IP3dSig() const;
  float DeltaEtaSeedCluster() const;
  float DeltaPhiSeedCluster() const;

  float PtError() const;

  float FBrem() const;

  float MvaID_Old() const;
  float MvaID_HZZ() const;
  float MvaID() const;
  float EnergyRegression() const;
  float EnergyRegressionErr() const;

  bool  PassConversionVeto() const;
  short ConversionMissHits() const;

  float ConversionDcot() const;
  float ConversionDist() const;
  float ConversionRadius() const;

  int TrackerLayersWithMeasurement() const;
  int NumberOfValidHits() const;
  int NumberOfValidPixelHits() const;
  int NumberOfValidTrackerHits() const;
  int NumberOfLostPixelHits() const;
  int NumberOfLostTrackerHits() const;


  TLorentzVector RegressionMomCombP4() const;

  float EffArea() const;

  //--------------------------
  // "set" methods ---------
  //--------------------------
  void AddTrack(TCElectron::Track);

  void SetNormalizedChi2Kf(float);

  void SetInverseEnergyMomentumDiff(float);

  void SetIP3d(float);
  void SetIP3dSig(float);
  void SetDeltaEtaSeedCluster(float);
  void SetDeltaPhiSeedCluster(float);

  void SetEoP(float);
  void SetEoPout(float);
  void SetESeedOverP(float);

  void SetPtError(float);

  void SetFBrem(float);

  void SetPassConversionVeto(bool);
  void SetConversionMissHits(short);

  void SetConversionDcot(float);
  void SetConversionDist(float);
  void SetConversionRadius(float);

  void SetTrackerLayersWithMeasurement(int);
  void SetNumberOfValidHits(int);
  void SetNumberOfValidPixelHits(int);
  void SetNumberOfValidTrackerHits(int);
  void SetNumberOfLostPixelHits(int);
  void SetNumberOfLostTrackerHits(int);


  void SetMvaID_Old(float);
  void SetMvaID_HZZ(float);
  void SetEnergyRegression(float);
  void SetEnergyRegressionErr(float);

  void SetRegressionMomCombP4(TLorentzVector tmpP4);

  void SetEffArea(float);

  ClassDef(TCElectron, 1);
};

#endif	/* _TCELECTRON_H */



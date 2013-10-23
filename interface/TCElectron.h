#ifndef _TCELECTRON_H
#define	_TCELECTRON_H

#include "TObject.h"
#include "TLorentzVector.h"
#include "TCEGamma.h"

class TCElectron : public TCEGamma {
 private:
  
  float _normChi2gsf;
  float _normChi2kf;

  float _ptError;
  
  float _fBrem;
  float _inverseEnergyMomentumDiff;

  float _mvaID;
  float _regEne;
  float _regErr;
  
  bool  _convVeto;
  short _convMissHits;
  
  
  int _trackerLayersWithMeasurement;
  int _numberOfValidHits;
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
  float _effArea;

 public:
  TCElectron();
  virtual ~TCElectron();
  
  // "get" methods -----------
  float NormalizedChi2() const;
  float NormalizedChi2Gsf() const;
  float NormalizedChi2Kf() const;

  float InverseEnergyMomentumDiff() const;
    
  float PtError() const;
  
  float FBrem() const;
  
  float MvaID() const; 
  float EnergyRegression() const; 
  float EnergyRegressionErr() const; 
  
  bool  ConversionVeto() const;
  short ConversionMissHits() const;
  
  
  int TrackerLayersWithMeasurement() const;
  int NumberOfValidHits() const;
  int NumberOfValidPixelHits() const;
  int NumberOfValidTrackerHits() const;
  int NumberOfLostPixelHits() const;
  int NumberOfLostTrackerHits() const;
  
  int  CutLevel(int lvl) const;
  bool PassID(int lvl) const;
  bool PassConversion(int lvl) const;
  bool PassIsolation(int lvl) const;
  
  TLorentzVector RegressionMomCombP4() const;

  float EffArea() const;

  //--------------------------
  // "set" methods ---------
  //--------------------------
  void SetNormalizedChi2Gsf(float);
  void SetNormalizedChi2Kf(float);
  
  void SetInverseEnergyMomentumDiff(float);

  void SetPtError(float);
  
  void SetFBrem(float);
  
  void SetConversionVeto(bool);
  void SetConversionMissHits(short);
  
  void SetTrackerLayersWithMeasurement(int);
  void SetNumberOfValidHits(int);
  void SetNumberOfValidPixelHits(int);
  void SetNumberOfValidTrackerHits(int);
  void SetNumberOfLostPixelHits(int);
  void SetNumberOfLostTrackerHits(int);
  
  
  void SetMvaID(float);
  void SetEnergyRegression(float);
  void SetEnergyRegressionErr(float);
  
  void SetCutLevel(int cut, int lvl);
  
  void SetRegressionMomCombP4(TLorentzVector tmpP4);
  
  void SetEffArea(float);

  ClassDef(TCElectron, 1);
};

#endif	/* _TCELECTRON_H */



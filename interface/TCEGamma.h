#ifndef _TCEGAMMA_H
#define	_TCEGAMMA_H

#include "TObject.h"
#include "TLorentzVector.h"
#include "TCPhysObject.h"

class TCEGamma : public TCPhysObject {


 public:
  struct CrystalInfo{
    int rawId;
    int ieta;
    int iphi;
    //int ix;
    //int iy;
    double energy;
    double time;
    double timeErr;
    int recoFlag;
  };


 private:
  
  bool  _isEB;                  // true if particle is in ECAL Barrel
  bool  _isEE;                  // true if particle is in ECAL Endcaps
  bool  _isInGap;
  
  float _hadOverEm;
  float _dPhiSuperCluster;
  float _dEtaSuperCluster;
  float _sigmaIetaIeta;
  float _eOverP;
  float _fBrem;
  float _r9;
  //float _mvaID;
  //float _regEne;
  //float _regErr;


  float _scEta;
  float _scPhi;

  // Notice that in case of electrons, the isolation variables are saved for the cone 0.4
  // while for the photons the cone is 0.3
  float _pfIsoCharged;
  float _pfIsoNeutral;
  float _pfIsoPhoton;

  //bool  _convVeto;
  //short _convMissHits;
  
  // crystal stuff
  vector<TCEGamma::CrystalInfo> _crysVect;
  int  _nCrystals; 

  
 public:
  TCEGamma();
  virtual ~TCEGamma();
  
   // "get" methods -----------
  vector<TCEGamma::CrystalInfo> GetCrystalVect() const;

  int   GetNCrystals() const;

  
  float HadOverEm() const;
  float DphiSuperCluster() const;
  float DetaSuperCluster() const;
  float SigmaIEtaIEta() const;
  float FBrem() const;
  float EOverP() const;
  
  float SCEta() const;
  float R9() const; 
  //float MvaID() const; 
  //float EnergyRegression() const; 
  //float EnergyRegressionErr() const; 
  
  //bool  ConversionVeto() const;
  //short ConversionMissHits() const;
  
  bool IsEB() const;
  bool IsEE() const;
  bool IsInGap() const;
  
  //bool PassConversion(int lvl) const;
  //bool PassIsolation(int lvl) const;
  
  //--------------------------
  // "set" methods ---------
  //--------------------------
  

  void AddCrystal(TCEGamma::CrystalInfo);
  void SetNCrystals(int);


  void SetHadOverEm(float h);
  void SetDphiSuperCluster(float dp);
  void SetDetaSuperCluster(float de);
  void SetSigmaIEtaIEta(float sieie);
  void SetEOverP(float e);
  void SetFBrem(float fb);
  
  
  
  void SetSCEta(float);
  
  //void SetConversionVeto(bool);
  //void SetConversionMissHits(short);
  
  void SetR9(float r);
  //void SetMvaID(float m);
  //void SetEnergyRegression(float e);
  //void SetEnergyRegressionErr(float e);
  
  void SetIsEB(bool b);
  void SetIsEE(bool b);
  void SetIsInGap(bool b);
  
  
  ClassDef(TCEGamma, 1);
};

#endif	/* _TCEGAMMA_H */



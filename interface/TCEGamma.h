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
  float _r9;
  
  //Superclaster shape variables. These should be coommon between electrons and photons.
  float _scEta;
  float _scPhi;
  float _scDeltaPhi;
  float _scDeltaEta;
  float _scSigmaIetaIeta;
  float _scSigmaIphiIphi;
  float _scEtaWidth;
  float _scPhiWidth;

  float _scEnergy;
  float _preShowerOverRaw;

  float _e1x5;
  float _e2x5;
  float _e5x5;
  

  //float _mvaID;
  //float _regEne;
  //float _regErr;


  // Notice that in case of electrons, the isolation variables are saved for the cone 0.4
  // while for the photons the cone is 0.3
  float _pfIsoCharged;
  float _pfIsoNeutral;
  float _pfIsoPhoton;

  //bool  _convVeto;
  
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
  
  float R9() const; 

  float SCEta() const;
  float SCPhi() const;
  float SCDeltaEta() const;
  float SCDeltaPhi() const;

  float SigmaIEtaIEta() const;
  float SigmaIPhiIPhi() const;

  float SCEtaWidth() const;
  float SCPhiWidth() const;

  float SCEnergy() const;
  float PreShowerOverRaw() const;
  float E1x5() const;
  float E2x5() const;
  float E5x5() const;

  float PfIsoCharged() const;
  float PfIsoNeutral() const;
  float PfIsoPhoton() const;

  //float MvaID() const; 
  //float EnergyRegression() const; 
  //float EnergyRegressionErr() const; 
  
  //bool  ConversionVeto() const;
  
  bool IsEB() const;
  bool IsEE() const;
  bool IsInGap() const;
  
  //bool PassConversion(int lvl) const;
  
  //--------------------------
  // "set" methods ---------
  //--------------------------
  
  void AddCrystal(TCEGamma::CrystalInfo);
  void SetNCrystals(int);

  void SetHadOverEm(float);

  void SetSCEta(float);
  void SetSCPhi(float);
  
  void SetSCDeltaEta(float);
  void SetSCDeltaPhi(float);

  void SetSigmaIEtaIEta(float);
  void SetSigmaIPhiIPhi(float);

  void SetSCEtaWidth(float);
  void SetSCPhiWidth(float);

  void SetSCEnergy(float);
  void SetPreShowerOverRaw(float);
  void SetE1x5(float);
  void SetE2x5(float);
  void SetE5x5(float);
  
  
  //void SetConversionVeto(bool);
  
  void SetR9(float);
  //void SetMvaID(float m);
  //void SetEnergyRegression(float e);
  //void SetEnergyRegressionErr(float e);
  
  void SetIsEB(bool);
  void SetIsEE(bool);
  void SetIsInGap(bool);
  
  void SetPfIsoCharged(float);
  void SetPfIsoNeutral(float);
  void SetPfIsoPhoton(float);
  
  ClassDef(TCEGamma, 1);
};

#endif	/* _TCEGAMMA_H */



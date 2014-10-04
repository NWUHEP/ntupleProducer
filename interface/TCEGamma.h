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
    CrystalInfo():
      rawId(-999),
      ieta(-999),
      iphi(-999),
      energy(-999),
      time(-999),
      timeErr(-999),
      recoFlag(-999)
    {}
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
  float _scSigmaIetaIphi;
  float _scSigmaIphiIphi;
  float _scEtaWidth;
  float _scPhiWidth;

  float _scRawEnergy;
  float _scEnergy;
  float _scPSEnergy;
  float _preShowerOverRaw;

  float _e1x3;
  float _e1x5;
  float _e2x2;
  float _e2x5;
  float _e2x5Max;
  float _e5x5;

  vector<float> _esEffSigmaRR;

  //float _mvaID;
  float _regEne;
  float _regErr;

  // Notice that in case of electrons, the isolation variables are saved for the cone 0.4
  // while for the photons the cone is 0.3
  float _pfIsoCharged;
  float _pfIsoNeutral;
  float _pfIsoPhoton;

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
  float SigmaIEtaIPhi() const;
  float SigmaIPhiIPhi() const;

  float SCEtaWidth() const;
  float SCPhiWidth() const;

  float SCRawEnergy() const;
  float SCEnergy() const;
  float SCPSEnergy() const;
  float PreShowerOverRaw() const;
  float E1x3() const;
  float E1x5() const;
  float E2x2() const;
  float E2x5() const;
  float E2x5Max() const;
  float E5x5() const;
  float E2OverE5() const;

  vector<float> ESEffSigmaRR() const;

  float PfIsoCharged() const;
  float PfIsoNeutral() const;
  float PfIsoPhoton() const;


  //float MvaID() const;
  float EnergyRegression() const;
  float EnergyRegressionErr() const;

  bool IsEB() const;
  bool IsEE() const;
  bool IsInGap() const;

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
  void SetSigmaIEtaIPhi(float);
  void SetSigmaIPhiIPhi(float);

  void SetSCEtaWidth(float);
  void SetSCPhiWidth(float);

  void SetSCRawEnergy(float);
  void SetSCEnergy(float);
  void SetSCPSEnergy(float);
  void SetPreShowerOverRaw(float);
  void SetE1x3(float);
  void SetE1x5(float);
  void SetE2x2(float);
  void SetE2x5(float);
  void SetE2x5Max(float);
  void SetE5x5(float);

  void SetESEffSigmaRR(float, float, float);

  void SetR9(float);
  //void SetMvaID(float m);
  void SetEnergyRegression(float);
  void SetEnergyRegressionErr(float);

  void SetIsEB(bool);
  void SetIsEE(bool);
  void SetIsInGap(bool);

  void SetPfIsoCharged(float);
  void SetPfIsoNeutral(float);
  void SetPfIsoPhoton(float);

  // print method
  virtual ostream& TCprint(ostream& out) const;

  ClassDef(TCEGamma, 1);
};

#endif	/* _TCEGAMMA_H */

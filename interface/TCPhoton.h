#ifndef _TCPHOTON_H
#define	_TCPHOTON_H

#include <memory>
#include "TObject.h"
#include "TLorentzVector.h"
#include "TArrayF.h"
#include "TCEGamma.h"
#include <vector>

class TCPhoton : public TCEGamma {
 public:
  struct FootprintRemoval{
    float chargediso;
    float chargediso_primvtx;
    float neutraliso;
    float photoniso;
    float chargediso_rcone;
    float chargediso_primvtx_rcone;
    float neutraliso_rcone;
    float photoniso_rcone;
    float eta_rcone;
    float phi_rcone;
    bool rcone_isOK;
    FootprintRemoval():
      chargediso(-1),
      chargediso_primvtx(-1),
      neutraliso(-1),
      photoniso(-1),
      chargediso_rcone(-1),
      chargediso_primvtx_rcone(-1),
      neutraliso_rcone(-1),
      photoniso_rcone(-1),
      eta_rcone(-1),
      phi_rcone(-1),
      rcone_isOK(false)
    {}
  };

private:
  // ID variables
  bool _trackVeto;
  bool _convVeto;
  int _nTrkSolidConeDR03;

  TCPhoton::FootprintRemoval _SCFootprintRemoval;

  // CiC iso vectors
  vector<float> _phoCiCPF4chgpfIso02;
  vector<float> _phoCiCPF4chgpfIso03;
  vector<float> _phoCiCPF4chgpfIso04;

 public:
    TCPhoton();
    virtual ~TCPhoton();

    // "get" methods -----------

    TCPhoton::FootprintRemoval GetSCFootprintRemovalStruct() const;

    bool  TrackVeto() const;
    bool  ConversionVeto() const;
    int NTrkSolidConeDR03() const;

    vector<float> CiCPF4chgpfIso02() const;
    vector<float> CiCPF4chgpfIso03() const;
    vector<float> CiCPF4chgpfIso04() const;

    // "set" methods ---------
    void SetSCFootprintRemovalStruct(TCPhoton::FootprintRemoval);

    void SetTrackVeto(bool);
    void SetConversionVeto(bool);
    void SetNTrkSolidConeDR03(int);

    void SetCiCPF4chgpfIso02(vector<float>);
    void SetCiCPF4chgpfIso03(vector<float>);
    void SetCiCPF4chgpfIso04(vector<float>);

    // print method
    virtual ostream& TCprint(ostream& out) const;

    ClassDef(TCPhoton, 1);
};

#endif



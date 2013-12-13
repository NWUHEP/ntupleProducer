#ifndef _TCPHOTON_H
#define	_TCPHOTON_H

#include <memory>
#include "TObject.h"
#include "TLorentzVector.h"
#include "TArrayF.h"
#include "TCEGamma.h"
#include <vector>

class TCPhoton : public TCEGamma {
private:

  // ID variables
  bool _trackVeto;
  bool _convVeto;
  int _nTrkSolidConeDR03;

 public:
    TCPhoton();
    virtual ~TCPhoton();

    // "get" methods -----------

    bool  TrackVeto() const;
    bool  ConversionVeto() const;
    int NTrkSolidConeDR03() const;

    // "set" methods ---------
    void SetTrackVeto(bool);
    void SetConversionVeto(bool);
    void SetNTrkSolidConeDR03(int);

    ClassDef(TCPhoton, 1);
};

#endif



#ifndef _TCPHOTON_H
#define	_TCPHOTON_H

#include <memory>
#include "TObject.h"
#include "TLorentzVector.h"
#include "TArrayF.h"
#include "TCEGamma.h"
#include <vector>

using namespace std;

class TCPhoton : public TCEGamma {
private:

  // ID variables
  //float _e2OverE9;
  bool  _trackVeto;
  
  bool    _convVeto;
  int _nTrkSolidConeDR03;

 public:
    TCPhoton();
    virtual ~TCPhoton();

    // "get" methods -----------

    //float E2OverE9() const; 
    bool  TrackVeto() const;
    bool  ConversionVeto() const;
    int NTrkSolidConeDR03() const;

    // "set" methods ---------
    //void SetE2OverE9(float);
    void SetTrackVeto(bool);
    void SetConversionVeto(bool);
    void SetNTrkSolidConeDR03(int);

    ClassDef(TCPhoton, 1);
};

#endif



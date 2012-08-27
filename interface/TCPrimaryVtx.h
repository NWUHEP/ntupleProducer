#ifndef _TCPRIMARYVTX_H
#define	_TCPRIMARYVTX_H

#include "TObject.h"
#include "TVector3.h"

class TCPrimaryVtx : public TVector3 {
private:

    float _nDof;
    float _chi2;
    bool  _isFake;
    int   _nTracks;
    float _sumPt2Trks;


public:
    TCPrimaryVtx();
    virtual ~TCPrimaryVtx();

    float NDof() const;
    float Chi2() const;
    bool  IsFake() const;
    int   Ntracks() const;
    float SumPt2Trks() const;

    // set methods
    void SetNDof(float n);
    void SetChi2(float chi2);
    void SetIsFake(bool isF);
    void SetNtracks(int nTrk);
    void SetSumPt2Trks(float sumPt2); 


    ClassDef(TCPrimaryVtx, 1);
};

#endif	/* _TCPRIMARYVTX_H */


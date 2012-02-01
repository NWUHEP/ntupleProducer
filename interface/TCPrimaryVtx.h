/* 
 * File:   TCPrimaryVtx.h
 * Author: Anton A.
 *
 * Created on May 21, 2010, 11:16 AM
 */

#ifndef _TCPRIMARYVTX_H
#define	_TCPRIMARYVTX_H

#include "TObject.h"
#include "TVector3.h"

class TCPrimaryVtx : public TObject {
private:

    TVector3 _position;
    float _nDof;
    float _chi2;
    bool  _isFake;
    int   _nTracks;
    float _sumPt2Trks;


public:
    TCPrimaryVtx();
    virtual ~TCPrimaryVtx();

    TVector3 Position() const;
    float NDof() const;
    float Chi2() const;
    bool  IsFake() const;
    int   Ntracks() const;
    float SumPt2Trks() const;

    // set methods
    void SetPosition(float x, float y, float z);
    void SetNDof(float n);
    void SetChi2(float chi2);
    void SetIsFake(bool isF);
    void SetNtracks(int nTrk);
    void SetSumPt2Trks(float sumPt2); 


    ClassDef(TCPrimaryVtx, 1);
};

#endif	/* _TCPRIMARYVTX_H */

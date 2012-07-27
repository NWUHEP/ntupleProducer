#ifndef _TCPHYSOBJECT_H
#define	_TCPHYSOBJECT_H

#include "TObject.h"
#include "TLorentzVector.h"
#include "TVector2.h"
#include "TVector3.h"
#include <map>
#include <utility>
#include <string>
#include <iostream>

using namespace std;

class TCPhysObject : public TLorentzVector {
    private:
        TVector3 _vtx;
        map<string, float> _IdMap;
        map<string, float> _IsoMap;
        int _charge;

    public:
        TCPhysObject();
        virtual ~TCPhysObject();

        // "get" methods -----------

        float IdMap(string key);
        float IsoMap(string key);
        TVector2 P2() const;
        TVector3 Vtx() const;
        int Charge() const;  
        float Dxy(TVector3 *primVtx) const;
        float Dz(TVector3 *primVtx) const;

        // "set" methods ---------
        void SetIdMap(string s, float v);
        void SetIsoMap(string s, float v);
        void SetVtx(float vx, float vy, float vz);
        void SetCharge(int c);  

        ClassDef(TCPhysObject, 1);
};

#endif	/* _TCPHYSOBJECT_H */

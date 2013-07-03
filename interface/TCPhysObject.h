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
        string _type;
        bool _isPF;
        bool _isReco;

    public:
        TCPhysObject();
        TCPhysObject(TLorentzVector p4, int charge);
        TCPhysObject(TLorentzVector p4, int charge, string type);
        virtual ~TCPhysObject();

        // "get" methods -----------

        float IdMap(string key);
        float IsoMap(string key);
        TVector2 P2() const;
        TVector3 Vtx() const;
        int Charge() const;  
        string Type() const;
        bool IsPF() const;
        bool IsReco() const;

        float Dxy(TVector3 *primVtx) const;
        float Dz(TVector3 *primVtx) const;

        // "set" methods ---------
        void SetP4(TLorentzVector p4);
        void SetIdMap(string s, float v);
        void SetIsoMap(string s, float v);
        void SetVtx(float vx, float vy, float vz);
        void SetCharge(int c);  
        void SetType(string s);
        void SetReco(bool);
        void SetPF(bool);

        ClassDef(TCPhysObject, 1);
};

#endif	/* _TCPHYSOBJECT_H */

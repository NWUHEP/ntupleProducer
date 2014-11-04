#ifndef _TCPHYSOBJECT_H
#define	_TCPHYSOBJECT_H


#include "TObject.h"
#include "TLorentzVector.h"
#include "TVector2.h"
#include "TVector3.h"
#include <map>
#include <algorithm>
#include <utility>
#include <string>
#include <iostream>
#include <assert.h>

using namespace std;

class TCPhysObject : public TLorentzVector {
    private:
        TVector3 _vtx;
        map<string, float> _IdMap;
        map<string, vector<string> > _Triggers;
        int  _charge;
        bool _isPF;
        bool _isTriggered;

    public:
        TCPhysObject();
        TCPhysObject(TLorentzVector p4, int charge);
        virtual ~TCPhysObject();

        // "get" methods -----------

        float IdMap(string key) const;
        TVector2 P2() const;
        TVector3 Vtx() const;
        int Charge() const;  
        bool IsPF() const;
        bool IsTriggered() const;

        float Dxy(TVector3 *primVtx) const;
        float Dz(TVector3 *primVtx) const;

        map<string, vector<string> > GetTriggers() const;


        // "set" methods ---------
        void SetP4(TLorentzVector p4);
        void SetIdMap(string s, float v);
        //void SetIsoMap(string s, float v);
        void SetVtx(float vx, float vy, float vz);
        void SetCharge(int c);  
        void SetPF(bool);
        void SetTriggered(bool);
        void AddTrigger(string hlt, string l3, const vector<string> hlts);

        // print method
        virtual ostream& TCprint(ostream& out) const;

        ClassDef(TCPhysObject, 1);
};

inline ostream& operator<<(ostream& os, const TCPhysObject& ph){ 
  return ph.TCprint(os);
}

#endif	/* _TCPHYSOBJECT_H */

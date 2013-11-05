#include "../interface/TCPhysObject.h"
#include "TCPhysObjectLinkDef.h"
//

TCPhysObject::TCPhysObject() {
    _isPF = _isReco = false;
}

TCPhysObject::TCPhysObject(TLorentzVector p4, int charge) {
    this->SetP4(p4);
    this->SetCharge(charge);

    _isPF = _isReco = false;
}

TCPhysObject::TCPhysObject(TLorentzVector p4, int charge, string type) {
    this->SetP4(p4);
    this->SetCharge(charge);
    this->SetType(type);

    _isPF = _isReco = false;
}

TCPhysObject::~TCPhysObject() {
}

// "get" methods -------------------------------------

using namespace std;

float TCPhysObject::IdMap(string key) const { 
    
    //Check that key is present in the id map
    try {
        string exception = "Can't find " + key + " in id map"; 
        if (_IdMap.count(key) == 0)
            throw exception;
    } catch (string ex) {
        cout << ex << endl;
    }

    return _IdMap.find(key)->second; 
}

float TCPhysObject::IsoMap(string key) const { 
    
    //Check that key is present in the iso map
    try {
        string exception = "Can't find " + key + " in isolation map"; 
        if (_IsoMap.count(key) == 0)
            throw exception;
    } catch (string ex) {
        cout << ex << endl;
    }

    return _IsoMap.find(key)->second; 
}

TVector2 TCPhysObject::P2() const {
    TVector2 v2(this->Px(), this->Py());
    return v2;
}

TVector3 TCPhysObject::Vtx() const { return _vtx; }
string TCPhysObject::Type() const { return _type; }
int TCPhysObject::Charge() const { return _charge; }
bool TCPhysObject::IsPF() const { return _isPF; }
bool TCPhysObject::IsReco() const { return _isReco; }

// "set" methods ---------------------------------------------

void TCPhysObject::SetP4(TLorentzVector p4) { this->SetPxPyPzE(p4.Px(), p4.Py(), p4.Pz(), p4.E()); } 
void TCPhysObject::SetIdMap(string s, float v){ _IdMap[s] = v; }
void TCPhysObject::SetIsoMap(string s, float v){ _IsoMap[s] = v; }

void TCPhysObject::SetVtx(float vx, float vy, float vz) {
    TVector3 v3(vx, vy, vz);
    _vtx = v3;
}

void TCPhysObject::SetCharge(int c){ _charge = c; }
void TCPhysObject::SetType(string s){ _type = s; }
void TCPhysObject::SetPF(bool p){ _isPF = p;}
void TCPhysObject::SetReco(bool r){ _isReco = r;}

// generally useful methods -----------------------------------

float TCPhysObject::Dxy(TVector3 *primVtx) const {
    //Calculating track dxy parameter 
    //wrt primary vertex d0 = - dxy                                                                                                                                               
    float vx = _vtx.X(), vy = _vtx.Y();
    float px = this->Px(), py = this->Py(), pt = this->Pt();
    float pvx = primVtx->X(), pvy = primVtx->Y();
    float ret =  (-(vx-pvx)*py + (vy-pvy)*px)/pt;
    return ret;
}

float TCPhysObject::Dz(TVector3 *primVtx) const {
    //Calculating track dz parameter wrt primary vertex                                                                                                        
    float vx = _vtx.X(), vy = _vtx.Y(), vz = _vtx.Z();
    float px = this->Px(), py = this->Py();
    float pz = this->Pz(), pt = this->Pt();
    float pvx = primVtx->X(), pvy = primVtx->Y(), pvz = primVtx->Z();
    float ret =  (vz-pvz)-((vx-pvx)*px +(vy-pvy)*py)/pt*(pz/pt);
    return ret;
}

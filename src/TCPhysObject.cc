#include "Higgs/ntupleProducer/interface/TCPhysObject.h"

TCPhysObject::TCPhysObject() {
}

TCPhysObject::~TCPhysObject() {
}

// "get" methods -------------------------------------

using namespace std;

float TCPhysObject::IdMap(string key) {
  return _IdMap[key];
}

float TCPhysObject::IsoMap(string key) {
  return _IsoMap[key];
}

TVector2 TCPhysObject::P2() const {
  TVector2 v2(this->Px(), this->Py());
  return v2;
}

TVector3 TCPhysObject::Vtx() const {
   return _vtx;
}


int TCPhysObject::Charge() const {
   return _charge;
}

float TCPhysObject::Dxy(TVector3 *primVtx) const {
  //Calculating track dxy parameter wrt primary vertex                                                                                                       
  //d0 = - dxy                                                                                                                                               
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

// "set" methods ---------------------------------------------

void TCPhysObject::SetIdMap(string s, float v){
  _IdMap[s] = v;
}

void TCPhysObject::SetIsoMap(string s, float v){
  _IsoMap[s] = v;
}

void TCPhysObject::SetVtx(float vx, float vy, float vz) {
   TVector3 v3(vx, vy, vz);
   _vtx = v3;
}

void TCPhysObject::SetCharge(int c){
   _charge = c;
}

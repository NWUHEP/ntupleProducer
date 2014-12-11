#include "../interface/TCTriggerObject.h"
#include "TCTriggerObjectLinkDef.h"

TCTriggerObject::TCTriggerObject() {
}

TCTriggerObject::~TCTriggerObject() {
}

void TCTriggerObject::SetId(int i) {
       _id = i;
}

void TCTriggerObject::SetHLTName(string s) {
       _HLTName = s;
}

void TCTriggerObject::SetModuleName(string s) {
       _moduleName = s;
}

void TCTriggerObject::SetEta(float e){
  _eta = e;
}

void TCTriggerObject::SetPhi(float p){
  _phi = p;
}

int TCTriggerObject::GetId() {
       return _id;
}

string TCTriggerObject::GetHLTName() {
       return _HLTName;
}

string TCTriggerObject::GetModuleName() {
       return _moduleName;
}

float TCTriggerObject::Eta() {
  return _eta;
}

float TCTriggerObject::Phi() {
  return _phi;
}

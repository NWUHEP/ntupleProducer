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

int TCTriggerObject::GetId() {
       return _id;
}

string TCTriggerObject::GetHLTName() {
       return _HLTName;
}

string TCTriggerObject::GetModuleName() {
       return _moduleName;
}

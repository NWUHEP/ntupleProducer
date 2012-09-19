#include "../interface/TCTriggerObject.h"

TCTriggerObject::TCTriggerObject() {
}

TCTriggerObject::~TCTriggerObject() {
}

void TCTriggerObject::SetId(int i) {
       _id = i;
}


int TCTriggerObject::GetId() {
       return _id;
}

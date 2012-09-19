#ifndef _TCTriggerObject_H
#define _TCTriggerObject_H

#include "TObject.h"
#include "TLorentzVector.h"

class TCTriggerObject : public TLorentzVector {
    private:
        int _id;

    public:
        TCTriggerObject();
        virtual ~TCTriggerObject();

        void SetId(int i);
        int GetId();

        ClassDef(TCTriggerObject, 1);
};

#endif

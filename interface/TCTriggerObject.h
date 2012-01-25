#ifndef _TCTriggerObject_H
#define _TCTriggerObject_H

#include "TObject.h"
#include "TLorentzVector.h"

class TCTriggerObject : public TObject {
    private:
        TLorentzVector p4;
        int id;

    public:
        TCTriggerObject();
        virtual ~TCTriggerObject();

        void setId(int i);
        void setP4(double px, double py, double pz, double energy);

        TLorentzVector P4();
        int getId();

        ClassDef(TCTriggerObject, 1);
};

#endif

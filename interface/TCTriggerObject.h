#ifndef _TCTriggerObject_H
#define _TCTriggerObject_H

#include "TObject.h"
#include "TLorentzVector.h"
#include <string>

using namespace std;

class TCTriggerObject : public TLorentzVector {
    private:
        int _id;
        string _HLTName;
        string _moduleName;

    public:
        TCTriggerObject();
        virtual ~TCTriggerObject();

        void SetId(int i);
        void SetHLTName(string s);
        void SetModuleName(string s);
        int GetId();
        string GetHLTName();
        string GetModuleName();

        ClassDef(TCTriggerObject, 1);
};

#endif

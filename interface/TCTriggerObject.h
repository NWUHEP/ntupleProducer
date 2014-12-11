#ifndef _TCTriggerObject_H
#define _TCTriggerObject_H

#include "TObject.h"
#include "TLorentzVector.h"
#include <string>

using namespace std;

class TCTriggerObject {
    private:
        int _id;
        string _HLTName;
        string _moduleName;
        float _eta;
        float _phi;
        
    public:
        TCTriggerObject();
        virtual ~TCTriggerObject();

        void SetId(int i);
        void SetHLTName(string s);
        void SetModuleName(string s);
        void SetEta(float e);
        void SetPhi(float p);
        int GetId();
        float Eta();
        float Phi();
        string GetHLTName();
        string GetModuleName();

        ClassDef(TCTriggerObject, 1);
};

#endif

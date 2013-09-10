#ifndef _TGENPARTICLE_H
#define _TGENPARTICLE_H

#include "TObject.h"
#include "TCPhysObject.h"
#include "TLorentzVector.h"
#include "TVector2.h"
#include "TVector3.h"
#include <vector>

class TCGenParticle : public TCPhysObject {
    private:
        TCGenParticle* mother;
        int PDGID;
        unsigned status;
        bool isParton_;


    public:
        TCGenParticle();
        virtual ~TCGenParticle();

        TCGenParticle* Mother();
        int GetPDGId();
        unsigned GetStatus();

        bool IsParton();

        void SetMother(TCGenParticle* m);
        void SetPDGId(int pdg_id);
        void SetStatus(unsigned s);

        void SetIsParton(bool a);

        ClassDef(TCGenParticle, 1);
};

#endif	

#ifndef _TGENPARTICLE_H
#define _TGENPARTICLE_H

#include "TCPhysObject.h"

class TCGenParticle : public TCPhysObject {
    private:
        TCGenParticle* _mother;
        int _PDGID;
        unsigned _status;
        bool _isParton;


    public:
        TCGenParticle();
        virtual ~TCGenParticle();

        TCGenParticle* Mother();
        //TCGenParticle* PrimaryAncestor();
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

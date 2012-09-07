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
        int mother;
        int grandmother;
        int PDGID;
        unsigned status;
        bool isParton_;

        //    std::vector<int> daughters;

    public:
        TCGenParticle();
        virtual ~TCGenParticle();

        int Mother();
        int Grandmother();
        int GetPDGId();
        unsigned GetStatus();

        bool IsParton();

        //    std::vector<int> GetDaughters();

        // "set" methods ---------
        //    void AddDaughter(int d);
        void SetMother(int m);
        void SetGrandmother(int g);
        void SetPDGId(int pdg_id);
        void SetStatus(unsigned s);

        void SetIsParton(bool a);

        ClassDef(TCGenParticle, 1);
};

#endif	

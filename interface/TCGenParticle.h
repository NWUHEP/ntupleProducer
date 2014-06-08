#ifndef _TGENPARTICLE_H
#define _TGENPARTICLE_H

#include "TCPhysObject.h"

class TCGenParticle : public TCPhysObject {
  private:
    int _momID;
    TCGenParticle* _mother;
    int _PDGID;
    unsigned _status;
    bool _isParton;


  public:
    TCGenParticle();
    virtual ~TCGenParticle();

    int MotherId() const;
    TCGenParticle* Mother();
    int GetPDGId() const;
    unsigned GetStatus() const;

    bool IsParton() const;

    void SetMotherId(int id);
    void SetMother(TCGenParticle* m);
    void SetPDGId(int pdg_id);
    void SetStatus(unsigned s);

    void SetIsParton(bool a);

    ClassDef(TCGenParticle, 2);

    // print method
    virtual ostream& TCprint(ostream& out) const;
};

#endif	

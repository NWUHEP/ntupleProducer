/* 
 * File:   TCJet.h
 * Author: Anton A.
 *
 * Created on April 30, 2010, 2:49 PM
 */

#ifndef _TCJET_H
#define	_TCJET_H

#include "TObject.h"
#include "TLorentzVector.h"
#include "TVector2.h"

class TCJet : public TObject {
private:
    TLorentzVector _p4;
    TVector3 _vtx;    
    //    TVector3 _assocPV;
    float _vtxSumPtFrac;
    float _vtxSumPt;
    float _vtxTrackFrac;
    float _vtxNTracks;
    unsigned int _vtxSumPtIndex;
    unsigned int _vtxCountIndex;

	float _jesUncertainty;
    float _jetCorr[8];
    bool  _jetCorrIsSet[8];

    float _chHadFrac;
    float _neuHadFrac;
    float _chEmFrac;
    float _neuEmFrac;

    unsigned int _numConstit;
    unsigned int _numChPart;

    // b tagging discriminators
    float _bDiscrTCHP;
    float _bDiscrTCHE;
    float _bDiscrSSVHE;
    float _bDiscrSSVHP;
    float _bDiscrJP;
    float _bDiscrJBP;
    float _bDiscrCSV;
    int   _jetFlavor;

public:
    TCJet();
    virtual ~TCJet();

    // "get" methods -----------

    TLorentzVector P4() const;
    TVector2 P2() const;
    float Et() const;
    float Pt() const;

    // accessors for corrected jets (the argument is the level of correction)
    // Note: in this implementation all lower-level corrections will be
    // applied as well. In the future, add overloaded methods for specifying individual 
    // corrections when/if needed
    TLorentzVector P4(unsigned int lvl) const;
    TVector2 P2(unsigned int lvl) const;
    float Et(unsigned int lvl) const;
    float Pt(unsigned int lvl) const;

    float TotalJetCorr(unsigned int lvl) const;

    float ChHadFrac() const;
    float NeuHadFrac() const;
    float ChEmFrac() const;
    float NeuEmFrac() const;

    unsigned int NumConstit() const;
    unsigned int NumChPart() const;

    TVector3 Vtx() const;
    float VtxSumPtFrac() const;
    float VtxSumPt() const;
    float VtxTrackFrac() const;
    int   VtxNTracks() const;
    unsigned int VtxSumPtIndex() const;
    unsigned int VtxCountIndex() const;
    //    TVector3 AssocVtx() const;
    bool  JetCorrIsSet(unsigned int lvl) const;
    float JetCorr(unsigned int lvl) const;
	float UncertaintyJES() const;

    // b tagging discriminators
    float BDiscrTCHP() const;
    float BDiscrTCHE() const;
    float BDiscrSSVHE() const;
    float BDiscrSSVHP() const;
    float BDiscrJP() const;
    float BDiscrJBP() const;
    float BDiscrCSV() const;
    int   JetFlavor() const;

    // "set" methods ---------
    void SetP4(TLorentzVector p4);
    void SetP4(float px, float py, float pz, float e);
    void SetVtx(float vx, float vy, float vz);    
    void SetVtxSumPtFrac(float f);
    void SetVtxSumPt(float p);
    void SetVtxTrackFrac(float f);
    void SetVtxNTracks(int n);
    void SetVtxSumPtIndex(unsigned int i);
    void SetVtxCountIndex(unsigned int i);

    void SetChHadFrac(float c);
    void SetNeuHadFrac(float n);
    void SetChEmFrac(float c);
    void SetNeuEmFrac(float n);
    void SetNumConstit(unsigned int n);
    void SetNumChPart(unsigned int n);
    void SetJetCorr(unsigned int lvl, float corr);
	void SetUncertaintyJES(float u);

    // b tagging discriminators
    // see the corresponding class members for description
    void SetBDiscrTCHE(float d);
    void SetBDiscrTCHP(float d);
    void SetBDiscrSSVHE(float d);
    void SetBDiscrSSVHP(float d);
    void SetBDiscrJP(float d);
    void SetBDiscrJBP(float d);
    void SetBDiscrCSV(float d);

    void SetJetFlavor(float f);

    ClassDef(TCJet, 1);

};

#endif	/* _TCJET_H */


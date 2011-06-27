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
    unsigned int _vtxIndex;

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
    //Track counting tag with N = 3: trackCountingHighPurBJetTags
    float _bDiscrTrkCountHiPure;
    //Track counting tag with N = 2: trackCountingHighEffBJetTags
    float _bDiscrTrkCountHiEff;
    //Simple secondary vertex b tag: simpleSecondaryVertexBJetTags
    float _bDiscrSecVtxSimple;
    //Combined SV b tag using likelihood ratios: combinedSVBJetTags
    float _bDiscrSecVtxL;
    //Combined SV b tag using MVA: combinedSVMVABJetTags
    float _bDiscrSecVtxMVA;

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
    unsigned int VtxIndex() const;
    //    TVector3 AssocVtx() const;
    bool  JetCorrIsSet(unsigned int lvl) const;
    float JetCorr(unsigned int lvl) const;
	 float UncertaintyJES() const;

    // b tagging discriminators
    float BDiscrTrkCountHiPure() const;
    float BDiscrTrkCountHiEff() const;
    float BDiscrSecVtxSimple() const;
    float BDiscrSecVtxL() const;
    float BDiscrSecVtxMVA() const;

    // "set" methods ---------
    void SetP4(TLorentzVector p4);
    void SetP4(float px, float py, float pz, float e);
    void SetVtx(float vx, float vy, float vz);    
    void SetVtxSumPtFrac(float vtxSumPtFrac);
    void SetVtxSumPt(float vtxSumPt);
    void SetVtxTrackFrac(float vtxTrackFrac);
    void SetVtxNTracks(int vtxNTracks);
    void SetVtxIndex(unsigned int vtxIndex);

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
    void SetBDiscrTrkCountHiPure(float d);
    void SetBDiscrTrkCountHiEff(float d);
    void SetBDiscrSecVtxSimple(float d);
    void SetBDiscrSecVtxL(float d);
    void SetBDiscrSecVtxMVA(float d);


    ClassDef(TCJet, 1);

};

#endif	/* _TCJET_H */


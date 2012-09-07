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
#include "TCPhysObject.h"

class TCJet : public TCPhysObject {
    private:
        float _vtxSumPtFrac;
        float _vtxSumPt;
        float _vtxTrackFrac;
        float _vtxNTracks;
        unsigned int _vtxSumPtIndex;
        unsigned int _vtxCountIndex;

        float _jesUncertainty;

        float _chHadFrac;
        float _neuHadFrac;
        float _chEmFrac;
        float _neuEmFrac;

        unsigned int _numConstit;
        unsigned int _numChPart;

        // b tagging discriminators
        map<string, float> _bDiscrMap;

        // Jet flavor
        int   _jetFlavor;

    public:
        TCJet();
        virtual ~TCJet();


        float ChHadFrac() const;
        float NeuHadFrac() const;
        float ChEmFrac() const;
        float NeuEmFrac() const;

        unsigned int NumConstit() const;
        unsigned int NumChPart() const;

        float VtxSumPtFrac() const;
        float VtxSumPt() const;
        float VtxTrackFrac() const;
        int   VtxNTracks() const;
        unsigned int VtxSumPtIndex() const;
        unsigned int VtxCountIndex() const;

        float UncertaintyJES() const;

        // b tagging discriminators
        float BDiscriminatorMap(string key);

        // Jet flavor
        int   JetFlavor() const;

        // "set" methods ---------
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
        void SetUncertaintyJES(float u);

        // b tagging discriminators
        // see the corresponding class members for description
        void SetBDiscriminatorMap(string key, float val);

        // Jet flavor
        void SetJetFlavor(float f);

        ClassDef(TCJet, 1);

};

#endif	/* _TCJET_H */


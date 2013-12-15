//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sat Dec 14 23:16:09 2013 by ROOT version 5.32/00
// from TTree eventTree/eventTree
// found on file: nuTuple.root
//////////////////////////////////////////////////////////

#ifndef myAnalyzer_h
#define myAnalyzer_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>

#include "../interface/TCPhysObject.h"
#include "../interface/TCJet.h"
#include "../interface/TCMET.h"
#include "../interface/TCEGamma.h"
#include "../interface/TCTrack.h"
#include "../interface/TCElectron.h"
#include "../interface/TCPhoton.h"
#include "../interface/TCMuon.h"
#include "../interface/TCTau.h"
#include "../interface/TCGenParticle.h"
#include "../interface/TCGenJet.h"
#include "../interface/TCPrimaryVtx.h"
#include "../interface/TCTriggerObject.h"

// Header file for the classes stored in the TTree if any.
#include <TClonesArray.h>
#include <TVector3.h>

// Fixed size dimensions of array or collections stored in the TTree if any.

class myAnalyzer : public TSelector {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain

   // Declaration of leaf types
   TClonesArray    *recoJets;
   TClonesArray    *recoElectrons;
   TClonesArray    *recoMuons;
   TClonesArray    *recoPhotons;
   TCMET           *recoMET;
   TCMET           *mva_MET;
   TCMET           *track_MET;
   TCMET           *T0MET;
   TCMET           *T2MET;
   TClonesArray    *genJets;
   TClonesArray    *genParticles;
   TClonesArray    *triggerObjects;
   TClonesArray    *primaryVtx;
   TVector3        *beamSpot;
   Int_t           nPUVertices;
   Float_t         nPUVerticesTrue;
   Bool_t          isRealData;
   UInt_t          runNumber;
   ULong64_t       eventNumber;
   UInt_t          lumiSection;
   UInt_t          bunchCross;
   Float_t         ptHat;
   Float_t         qScale;
   Float_t         evtWeight;
   Float_t         rhoFactor;
   Float_t         rho25Factor;
   Float_t         rhoMuFactor;
   ULong64_t       triggerStatus;
   UInt_t          hltPrescale[64];
   Bool_t          NoiseFilters_isScraping;
   Bool_t          NoiseFilters_isNoiseHcalHBHE;
   Bool_t          NoiseFilters_isNoiseHcalLaser;
   Bool_t          NoiseFilters_isNoiseEcalTP;
   Bool_t          NoiseFilters_isNoiseEcalBE;
   Bool_t          NoiseFilters_isCSCTightHalo;
   Bool_t          NoiseFilters_isCSCLooseHalo;
   Bool_t          NoiseFilters_isNoiseTracking;
   Bool_t          NoiseFilters_isNoiseEEBadSc;
   Bool_t          NoiseFilters_isNoisetrkPOG1;
   Bool_t          NoiseFilters_isNoisetrkPOG2;
   Bool_t          NoiseFilters_isNoisetrkPOG3;

   // List of branches
   TBranch        *b_recoJets;   //!
   TBranch        *b_recoElectrons;   //!
   TBranch        *b_recoMuons;   //!
   TBranch        *b_recoPhotons;   //!
   TBranch        *b_recoMET;   //!
   TBranch        *b_mva_MET;   //!
   TBranch        *b_track_MET;   //!
   TBranch        *b_T0MET;   //!
   TBranch        *b_T2MET;   //!
   TBranch        *b_genJets;   //!
   TBranch        *b_genParticles;   //!
   TBranch        *b_triggerObjects;   //!
   TBranch        *b_primaryVtx;   //!
   TBranch        *b_beamSpot;   //!
   TBranch        *b_nPUVertices;   //!
   TBranch        *b_nPUVerticesTrue;   //!
   TBranch        *b_isRealData;   //!
   TBranch        *b_runNumber;   //!
   TBranch        *b_eventNumber;   //!
   TBranch        *b_lumiSection;   //!
   TBranch        *b_bunchCross;   //!
   TBranch        *b_ptHat;   //!
   TBranch        *b_qScale;   //!
   TBranch        *b_evtWeight;   //!
   TBranch        *b_rhoFactor;   //!
   TBranch        *b_rho25Factor;   //!
   TBranch        *b_rhoMuFactor;   //!
   TBranch        *b_triggerStatus;   //!
   TBranch        *b_hltPrescale;   //!
   TBranch        *b_NoiseFilters;   //!

   myAnalyzer(TTree * /*tree*/ =0) : fChain(0) { }
   virtual ~myAnalyzer() { }
   virtual Int_t   Version() const { return 2; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();

   ClassDef(myAnalyzer,0);
};

#endif

#ifdef myAnalyzer_cxx
void myAnalyzer::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   recoJets = 0;
   recoElectrons = 0;
   recoMuons = 0;
   recoPhotons = 0;
   recoMET = 0;
   mva_MET = 0;
   track_MET = 0;
   T0MET = 0;
   T2MET = 0;
   genJets = 0;
   genParticles = 0;
   triggerObjects = 0;
   primaryVtx = 0;
   beamSpot = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("recoJets", &recoJets, &b_recoJets);
   fChain->SetBranchAddress("recoElectrons", &recoElectrons, &b_recoElectrons);
   fChain->SetBranchAddress("recoMuons", &recoMuons, &b_recoMuons);
   fChain->SetBranchAddress("recoPhotons", &recoPhotons, &b_recoPhotons);
   fChain->SetBranchAddress("recoMET", &recoMET, &b_recoMET);
   fChain->SetBranchAddress("mva_MET", &mva_MET, &b_mva_MET);
   fChain->SetBranchAddress("track_MET", &track_MET, &b_track_MET);
   fChain->SetBranchAddress("T0MET", &T0MET, &b_T0MET);
   fChain->SetBranchAddress("T2MET", &T2MET, &b_T2MET);
   fChain->SetBranchAddress("genJets", &genJets, &b_genJets);
   fChain->SetBranchAddress("genParticles", &genParticles, &b_genParticles);
   fChain->SetBranchAddress("triggerObjects", &triggerObjects, &b_triggerObjects);
   fChain->SetBranchAddress("primaryVtx", &primaryVtx, &b_primaryVtx);
   fChain->SetBranchAddress("beamSpot", &beamSpot, &b_beamSpot);
   fChain->SetBranchAddress("nPUVertices", &nPUVertices, &b_nPUVertices);
   fChain->SetBranchAddress("nPUVerticesTrue", &nPUVerticesTrue, &b_nPUVerticesTrue);
   fChain->SetBranchAddress("isRealData", &isRealData, &b_isRealData);
   fChain->SetBranchAddress("runNumber", &runNumber, &b_runNumber);
   fChain->SetBranchAddress("eventNumber", &eventNumber, &b_eventNumber);
   fChain->SetBranchAddress("lumiSection", &lumiSection, &b_lumiSection);
   fChain->SetBranchAddress("bunchCross", &bunchCross, &b_bunchCross);
   fChain->SetBranchAddress("ptHat", &ptHat, &b_ptHat);
   fChain->SetBranchAddress("qScale", &qScale, &b_qScale);
   fChain->SetBranchAddress("evtWeight", &evtWeight, &b_evtWeight);
   fChain->SetBranchAddress("rhoFactor", &rhoFactor, &b_rhoFactor);
   fChain->SetBranchAddress("rho25Factor", &rho25Factor, &b_rho25Factor);
   fChain->SetBranchAddress("rhoMuFactor", &rhoMuFactor, &b_rhoMuFactor);
   fChain->SetBranchAddress("triggerStatus", &triggerStatus, &b_triggerStatus);
   fChain->SetBranchAddress("hltPrescale", hltPrescale, &b_hltPrescale);
   fChain->SetBranchAddress("NoiseFilters", &NoiseFilters_isScraping, &b_NoiseFilters);
}

Bool_t myAnalyzer::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

#endif // #ifdef myAnalyzer_cxx

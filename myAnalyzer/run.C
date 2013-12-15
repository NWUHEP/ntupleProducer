void run(){
  //Let root know where do we keep our .h and .cc files:
  gROOT->SetMacroPath(".:../src/:../interface/");
  
  //Need to load the libraries for all container classes
  gROOT->LoadMacro("TCPhysObject.cc+");
  gROOT->LoadMacro("TCJet.cc+");
  gROOT->LoadMacro("TCMET.cc+");
  gROOT->LoadMacro("TCEGamma.cc+");
  gROOT->LoadMacro("TCTrack.cc+");
  gROOT->LoadMacro("TCElectron.cc+");
  gROOT->LoadMacro("TCPhoton.cc+");
  gROOT->LoadMacro("TCMuon.cc+");
  gROOT->LoadMacro("TCTau.cc+");
  gROOT->LoadMacro("TCGenJet.cc+");
  gROOT->LoadMacro("TCGenParticle.cc+");
  gROOT->LoadMacro("TCPrimaryVtx.cc+");
  gROOT->LoadMacro("TCTriggerObject.cc+");
  
  TChain* fChain = new TChain("ntupleProducer/eventTree");
  
  fChain->Add("nuTuple.root"); 
  //fChain->Add("nuTuple_2.root"); 
  //etc. you can add loop here to add all the files you want to run over
  
  
  //And here is how we finally run it:
  fChain->Process("myAnalyzer.C+", "options");
}

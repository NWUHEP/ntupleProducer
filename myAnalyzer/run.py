#!/usr/bin/env python
from ROOT import *
gROOT.SetMacroPath(".:../src/:../interface/");

gROOT.LoadMacro("TCPhysObject.cc+");
gROOT.LoadMacro("TCJet.cc+");
gROOT.LoadMacro("TCMET.cc+");
gROOT.LoadMacro("TCEGamma.cc+");
gROOT.LoadMacro("TCTrack.cc+");
gROOT.LoadMacro("TCElectron.cc+");
gROOT.LoadMacro("TCPhoton.cc+");
gROOT.LoadMacro("TCMuon.cc+");
gROOT.LoadMacro("TCTau.cc+");
gROOT.LoadMacro("TCGenJet.cc+");
gROOT.LoadMacro("TCGenParticle.cc+");
gROOT.LoadMacro("TCPrimaryVtx.cc+");
gROOT.LoadMacro("TCTriggerObject.cc+");

fChain = TChain("ntupleProducer/eventTree");

fChain.Add("nuTuple.root"); 

fChain.Process("myAnalyzer.C+", "options");

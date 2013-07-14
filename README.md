The Ntuple Producer code
========================
This is a CMSSW code for creating small ROOT ntuples from CMS data/MC samples.

Instructions for Users
---------------------
 * Set up the environment
```
  setenv SCRAM_ARCH slc5_amd64_gcc462
  cmsrel CMSSW_5_3_8
  cd CMSSW_5_3_8/src
  cmsenv
```

 * Check out the ntuple producer code
```git clone git@github.com:NWUHEP/ntupleProducer.git NWU/ntupleProducer```

 * Additional packages for running PF2PAT and MET corrections in 53X
```
  addpkg CommonTools/ParticleFlow   V00-03-16
  addpkg RecoParticleFlow/PFProducer V15-01-11 
  addpkg DataFormats/METReco V03-03-11-01 
  addpkg JetMETCorrections/Type1MET V04-06-09-02
  addpkg CommonTools/RecoUtils V00-01-04
  addpkg DataFormats/ParticleFlowCandidate V15-03-04-01
  addpkg DataFormats/TrackReco V10-02-02-01
  addpkg DataFormats/VertexReco V02-00-04-01
```

 * Do we need this?? Since we don't use PAT anymore, perhaps, those are not needed:
```
  addpkg PhysicsTools/PatAlgos V08-09-52  
  cvs up -r V08-09-07-05 PhysicsTools/PatAlgos/python/patTemplate_cfg.py   
  addpkg PhysicsTools/PatUtils V03-09-28

  addpkg DataFormats/StdDictionaries V00-02-14
  addpkg FWCore/GuiBrowsers V00-00-70
  addpkg RecoMET/METProducers V03-03-12-02
  addpkg DataFormats/PatCandidates V06-05-06-07
```            

 * Noise, MET, tracker filters, following instructions here: CMS/MissingETOptionalFilters:
```
  cvs co -r V00-00-13 RecoMET/METFilters
  cvs co -r V00-00-08 RecoMET/METAnalyzers
  cvs co -r V00-03-23 CommonTools/RecoAlgos
  cvs co -r V01-00-11-01 DPGAnalysis/Skims
  cvs co -r V00-11-17 DPGAnalysis/SiStripTools
  cvs co -r V00-00-08 DataFormats/TrackerCommon
  cvs co -r V01-09-05 RecoLocalTracker/SubCollectionProducers
```

 * For calculation PFIso for EGamma objects,
```
  cvs co -r V00-00-30-02 -d EGamma/EGammaAnalysisTools UserCode/EGamma/EGammaAnalysisTools
  cvs up -r 1.13 EGamma/EGammaAnalysisTools/interface/PFIsolationEstimator.h
  cvs up -r 1.20 EGamma/EGammaAnalysisTools/src/PFIsolationEstimator.cc
```

 * Also this for the electron MVA ID and Regression
```
  cvs co -r HCP2012_V05 EgammaAnalysis/ElectronTools
  cvs co -r V09-00-01   RecoEgamma/EgammaTools
```

 * Finally, compile this mess (takes a while... coffee time!)  
 ```scram b -j 9```


Check out this twiki for more details: [CMS/UserCodeNWUntupleProducer][1]


Instructions for Developers
--------------------------

### Tagging policy


[1]: https://twiki.cern.ch/twiki/bin/view/CMS/UserCodeNWUntupleProducer
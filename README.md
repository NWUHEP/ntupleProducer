The Ntuple Producer code
========================
This is a CMSSW code for creating small ROOT ntuples from CMS data/MC samples.
Check out this twiki for more details: [CMS/UserCodeNWUntupleProducer][1]

Instructions for Users
---------------------
 * Set up the environment
```
  setenv SCRAM_ARCH slc5_amd64_gcc462
  cmsrel CMSSW_5_3_8
  cd CMSSW_5_3_8/src
  cmsenv
```

 * Check out the ntuple producer code and then the specific tag of the code that is known to work
```
 git clone git@github.com:NWUHEP/ntupleProducer.git NWU/ntupleProducer <-- This does not seem to work for me...
 cd NWU/ntupleProducer
 git checkout v6.3
 cd ../..
```

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

 * Do we need this?? Since we do not use PAT anymore, perhaps, those are not needed:
```
  addpkg PhysicsTools/PatAlgos V08-09-52  
  cvs up -r V08-09-07-05 PhysicsTools/PatAlgos/python/patTemplate_cfg.py   
  addpkg PhysicsTools/PatUtils V03-09-28

  addpkg DataFormats/StdDictionaries V00-02-14
  addpkg FWCore/GuiBrowsers V00-00-70
  addpkg RecoMET/METProducers V03-03-12-02
  addpkg DataFormats/PatCandidates V06-05-06-07
```            

 * Noise, MET, tracker filters, following instructions here: [CMS/MissingETOptionalFilters][2]:
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
```
 scram b -j 9
```

Once compiled, we are ready to run it
### Runnning the code
```
  cd test
  cmsRun ntupleProducer_cfg.py
```

*NB* 
For now, the ntuples require that there be at least one lepton with pT > 10 GeV in order for an event to be saved. In the case that this is not desired (for instance, in jet or photon based studies), you should modify the following line in the ntupleProducer.cc file
```c++ 
if (eleCount > 0 || muCount > 0) eventTree -> Fill();
```

In addition to this, there is a flag in the configuration file, ntupleProducer_cfg.py

#### Running with CRAB


Instructions for Developers
--------------------------

 * First, make sure you are on master branch and have the latest code:
```
  git checkout master
  git pull
```

 * Then create a new branch and swich to it:
```
  git branch dev1
  git checkout dev1
```

 * Now you can make any changes you want. Once you are done, commit it and push your branch.
```
  git commit -a
  git push origin dev1
```

 * When you are satisfied with you new code, merge it with master branch. For that:
```
  git checkout master
  git merge dev1
  git push
```

If the changes do not conflict, you are done. 
If there are conflicts, markers will be left in the problematic files showing the conflict; `git diff` will show this. 
Once you have edited the files to resolve the conflicts, `git commit -a`.
 
### Tagging policy
At any time you can tag your code, and push your tags to remote:
```
  git tag -a test1 -m "my tag"
  git push origin --tags
```
You can use any tags you want, later those can be deleted.

For the global production though, we should stick with a tagging convention.
Tags should be **vX.Y** and I am starting them with **v6.1**. Such that the tag corresponds to the **nutuple_v6** name 
of ntuple production. 
If the new code significantly changes the format of the ntuples (substantial changes to class definitions etc.) then the first number of a tag should be incremented 
(to v7.1 etc.) and the the ntuple production path-name should changed correspondingly.  Otherwise, incremental changes should be reflected in changes to the second digit


[1]: https://twiki.cern.ch/twiki/bin/view/CMS/UserCodeNWUntupleProducer
[2]: https://twiki.cern.ch/twiki/bin/view/CMS/MissingETOptionalFilters

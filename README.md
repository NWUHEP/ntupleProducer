The Ntuple Producer code
========================
This is a CMSSW code for creating small ROOT ntuples from CMS data/MC samples.
Check out this twiki for more details: [CMS/UserCodeNWUntupleProducer][1]

Instructions for Users
---------------------
 * Set up the environment
```
  setenv SCRAM_ARCH slc5_amd64_gcc462
  cmsrel CMSSW_5_3_8_patch1
  cd CMSSW_5_3_8_patch1/src
  cmsenv
```

 * Check out the ntuple producer code and then the specific tag of the code that is known to work
```
 git clone https://github.com/NWUHEP/ntupleProducer NWU/ntupleProducer
 cd NWU/ntupleProducer
 git checkout v8.0
 cd ../..
```

 * Met recipes, according to [2] and [3]:
```
  git cms-addpkg PhysicsTools/PatAlgos
  git cms-merge-topic -u vadler:53X-tagset133511
  git cms-addpkg PhysicsTools/PatUtils
  git cms-merge-topic -u TaiSakuma:53X-met-130910-01
```
 
* Met filters according to [4]:
```
  cvs co -r V00-00-13-01 RecoMET/METFilters
  cvs co -r V00-00-08 RecoMET/METAnalyzers
  cvs co -r V00-03-23 CommonTools/RecoAlgos
  cvs co -r V01-00-11-01 DPGAnalysis/Skims
  cvs co -r V00-11-17 DPGAnalysis/SiStripTools
  cvs co -r V00-00-08 DataFormats/TrackerCommon
  cvs co -r V01-09-05 RecoLocalTracker/SubCollectionProducers
```

 * Egamma tools from [5]:
```
  cvs co -r V00-00-09 EgammaAnalysis/ElectronTools
  cvs co -r V09-00-01 RecoEgamma/EgammaTools
  cd EgammaAnalysis/ElectronTools/data
  cat download.url | xargs wget
  cd ../../../
```

 * Track MET Code [need a ref]:
```
  cvs co -r 1.2 RecoMET/METProducers/src/ParticleFlowForChargedMETProducer.cc
  cvs co -r 1.1 RecoMET/METProducers/src/TrackMETProducer.cc
  cvs co -r 1.2 RecoMET/METProducers/interface/ParticleFlowForChargedMETProducer.h
  cvs co -r 1.1 RecoMET/METProducers/interface/TrackMETProducer.h
  cvs up -r 1.17 RecoMET/METProducers/src/SealModule.cc
  cvs co -r 1.1 RecoMET/METProducers/python/TrackMET_cfi.py
  cvs co -r 1.2 RecoMET/METProducers/python/pfChargedMET_cfi.py
```

 * MVA MET Code [need a ref]:
```
  cvs co -r METPU_5_3_X_v12 JetMETCorrections/METPUSubtraction
  cvs co -r HEAD -d pharrisTmp UserCode/pharris/MVAMet/data
  cp  -d pharrisTmp/*June2013*.root           JetMETCorrections/METPUSubtraction/data/
  cp  -d pharrisTmp/*Dec2012*.root           JetMETCorrections/METPUSubtraction/data/
  rm -rf pharrisTmp
  cvs co -r METPU_5_3_X_v4 RecoJets/JetProducers
  cvs up -r HEAD RecoJets/JetProducers/data/
  cvs up -r HEAD RecoJets/JetProducers/python/PileupJetIDCutParams_cfi.py                     
  cvs up -r HEAD RecoJets/JetProducers/python/PileupJetIDParams_cfi.py                     
  cvs up -r HEAD RecoJets/JetProducers/python/PileupJetID_cfi.py     
  cvs co -r b5_3_X_cvMEtCorr_2013Feb22            DataFormats/METReco
  cvs co -r V05-00-16                             DataFormats/JetReco
  cvs co -r V01-04-25                             RecoTauTag/RecoTau 
  cvs co -r V03-04-07                             RecoMET/METAlgorithms
  cvs co -r V01-04-13                             RecoTauTag/Configuration
```

 * Files that needs to be updated:
```
  cvs co -r V03-04-07 DataFormats/METReco/interface/CorrMETData.h
  cvs co -r HEAD JetMETCorrections/Type1MET/plugins/Type0PFMETcorrInputProducer.h
  cvs co -r HEAD JetMETCorrections/Type1MET/plugins/Type0PFMETcorrInputProducer.cc
```

 * Patches to checked folders:
```
  cp NWU/ntupleProducer/patches/pfMETCorrections_cff.py JetMETCorrections/Type1MET/python/pfMETCorrections_cff.py
  cp NWU/ntupleProducer/patches/mvaPFMET_leptons_data_cff.py JetMETCorrections/METPUSubtraction/python/mvaPFMET_leptons_data_cff.py
  cp NWU/ntupleProducer/patches/mvaPFMET_leptons_cff.py JetMETCorrections/METPUSubtraction/python/mvaPFMET_leptons_cff.py
  cp NWU/ntupleProducer/patches/mvaPFMET_leptons_cfi.py JetMETCorrections/METPUSubtraction/python/mvaPFMET_leptons_cfi.py
  cp NWU/ntupleProducer/patches/PATMHTProducer.h PhysicsTools/PatAlgos/plugins/PATMHTProducer.h
  cp NWU/ntupleProducer/patches/classes.h DataFormats/StdDictionaries/src/classes.h
  cp NWU/ntupleProducer/patches/classes_def.xml DataFormats/StdDictionaries/src/classes_def.xml
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
it assumes you are running over an MC sample. If you want to run on data, do:
```
  cmsRun ntupleProducer_cfg.py isRealData=1
``` 
that will set up an appropriate global tag etc.

*NB* 
By defualt, the ntuples require that there be at least one muon(electron) with pT > 3(5) GeV in order for an event to be saved. 
In the case that this is not desired (for instance, in jet or photon based studies), 
you should modify switch off the ```skimLeptons``` option in ntupleProducer_cfg.py

In addition to this, there are various flags the configuration file, ntupleProducer_cfg.py, that allow to save/not save certain objects (muons, jets, etc). All are saved by default.  

#### Running with CRAB
Look into ```crabNtuples_MC.cfg``` and ```crabNtuples_Data.cfg``` scripts.

Will incorporate multicrab soon.

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
(to v7.1 etc.) and the ntuple production path-name should be changed correspondingly.  Otherwise, incremental changes should be reflected in changes to the second digit.


[1]: https://twiki.cern.ch/twiki/bin/view/CMS/UserCodeNWUntupleProducer
[3]: https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMETRecipe53X
[2]: https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMetAnalysis
[4]: https://twiki.cern.ch/twiki/bin/view/CMS/MissingETOptionalFilters
[5]: https://twiki.cern.ch/twiki/bin/view/CMS/MultivariateElectronIdentification

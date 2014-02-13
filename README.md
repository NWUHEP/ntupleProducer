The Ntuple Producer code
========================
This is a CMSSW code for creating small ROOT ntuples from CMS data/MC samples.
Check out this twiki for more details: [CMS/UserCodeNWUntupleProducer][1]

Instructions for Users
---------------------
 * Set up the environment
```
  setenv SCRAM_ARCH slc5_amd64_gcc462
  setenv CVSROOT :ext:<cern-user-account>@lxplus.cern.ch:/afs/cern.ch/user/c/cvscmssw/public/CMSSW
  cmsrel CMSSW_5_3_13_patch3
  cd CMSSW_5_3_13_patch3/src
  cmsenv
```
Replace the ```<cern-user-account>``` with your CERN account.
Since the new recipe for CVS connection is done through ssh to lxplus you will have to type your CERN password
every time when checkout from CVS. It is inconveniet, but we have to live with it for now. 

 * Met recipes, according to [workbook][2] and [met-recipe][3]:
```
  git cms-addpkg PhysicsTools/PatAlgos
  git cms-merge-topic 1472
  git cms-merge-topic -u TaiSakuma:53X-met-131120-01
  scram b -j 9
```

 * Met filters according to [MissingETOptionalFilters][4]:
```
  cvs co -r V00-00-13-01 RecoMET/METFilters
  cvs co -r V00-03-23 CommonTools/RecoAlgos
  cvs co -r V01-00-11-01 DPGAnalysis/Skims
  cvs co -r V00-11-17 DPGAnalysis/SiStripTools
  cvs co -r V00-00-08 DataFormats/TrackerCommon
  cvs co -r V01-09-05 RecoLocalTracker/SubCollectionProducers
  scram b -j 9
```

 * Ecal tools so we can get more photon variables for photon MVA:
```
  git cms-addpkg RecoEcal/EgammaCoreTools
```

 * Egamma tools from [MVA Electron ID][5]:
```
  cvs co -r V00-00-09 EgammaAnalysis/ElectronTools
  cvs co -r V09-00-01 RecoEgamma/EgammaTools
  cd EgammaAnalysis/ElectronTools/data
  cat download.url | xargs wget
  cd ../../../
  scram b -j 9
```

 * MVA MET Code (Just for PU Jet ID) [Jet PU ID][9]:
```
  cvs co -r METPU_5_3_X_v4 RecoJets/JetProducers
  cvs up -r HEAD RecoJets/JetProducers/data/
  cvs up -r HEAD RecoJets/JetProducers/python/PileupJetIDCutParams_cfi.py
  cvs up -r HEAD RecoJets/JetProducers/python/PileupJetIDParams_cfi.py
  cvs up -r HEAD RecoJets/JetProducers/python/PileupJetID_cfi.py
  cvs co -r V05-00-16 DataFormats/JetReco
  scram b -j 9
```

 * Extra code (for boosted Z->ee isolation), (-d option not working since new CVS directory) following [boostedZ][6] and [heep][7] twikies:
```
  mkdir TSWilliams
  mkdir SHarper
  cvs co -r V00-02-03 UserCode/TSWilliams/BstdZee/BstdZeeTools
  cvs co -r V00-09-03 UserCode/SHarper/HEEPAnalyzer
  mv UserCode/TSWilliams/BstdZee/BstdZeeTools TSWilliams/.
  mv UserCode/SHarper/HEEPAnalyzer SHarper/.
  rm -r UserCode
```

 * PF footprint removal [Supercluster footprint removal twiki][8]:
```
  git clone https://github.com/peruzzim/SCFootprintRemoval.git PFIsolation/SuperClusterFootprintRemoval
  cd  PFIsolation/SuperClusterFootprintRemoval
  git checkout V01-06
  cd ../..
  scram b -j 9
```

 * Now check out the ntuple producer code and then the specific tag/branch of the code that is known to work
```
 git clone https://github.com/NWUHEP/ntupleProducer NWU/ntupleProducer
 cd NWU/ntupleProducer
 git checkout v9.6
 cd ../..
 scram b -j 9
```

Once compiled, we are ready to run it
### Runnning the code
```
  cd NWU/ntupleProducer/test
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
you should switch off the ```skimLeptons``` option in ntupleProducer_cfg.py

In addition to this, there are various flags the configuration file, ntupleProducer_cfg.py, that allow to save/not save certain objects (muons, jets, etc). All are saved by default.

#### Running with CRAB
Look into ```crabNtuples_MC.cfg``` and ```crabNtuples_Data.cfg``` scripts.

Multicrab added, will update with instructions soon.

#### Checking Output
After CRAB claims that your jobs are finished with exit codes 0 0, you will want to double check because it lies and large jobs tend to
have a few extra or missing files.

Run the following command:
```
  ./find_goodfiles.py -c Path/To/CrabDir -q
```
This will check that all the jobs listed in the crab xml files are actually in your output area, and that your output area contains no
extra or duplicate files.  If it does, the script will tell you what needs to be rerun or what needs to be deleted.

Instructions for Developers
--------------------------

 * First, make sure you are on master branch and have the latest code:
```
  git checkout master
  git pull
```

 * Then create a new branch and swich to it:
```
  git branch dev-username
  git checkout dev-username
```

 * Now you can make any changes you want. Once you are done, commit it and push your branch.
```
  git commit -a
  git push origin dev-username
```

 * When you are satisfied with you new code, merge it with master branch. For that:
```
  git checkout master
  git merge dev-username
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
[2]: https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMetAnalysis
[3]: https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMETRecipe53X
[4]: https://twiki.cern.ch/twiki/bin/view/CMS/MissingETOptionalFilters
[5]: https://twiki.cern.ch/twiki/bin/view/CMS/MultivariateElectronIdentification
[6]: https://twiki.cern.ch/twiki/bin/viewauth/CMS/BoostedZToEEModIso
[7]: https://twiki.cern.ch/twiki/bin/view/CMS/HEEPSelector
[8]: https://twiki.cern.ch/twiki/bin/viewauth/CMS/SuperClusterFootprintRemoval
[9]: https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupJetID

Several packages need to be check out for the ntuple code to run properly.
This is true as of CMSSW_4_2_4.  

--- Ecal noise cleaning
cvs co -r V110523_BE -d DataFormats/AnomalousEcalDataFormats UserCode/csander/DataFormats/AnomalousEcalDataFormats
cvs co -r V110523_BE -d PhysicsTools/EcalAnomalousEventFilter UserCode/csander/PhysicsTools/EcalAnomalousEventFilter
cvs co -r V110523_BE -d Sandbox UserCode/csander/Sandbox
cvs co -r V19MAY2011_v3 JetMETAnalysis/ecalDeadCellTools
#--- MET corrections
cvs co -r V04-04-04 JetMETCorrections/Type1MET

 Log.           It's better than bad it's good.

2012 Sep 07
  * Tagged B11-01 for Nate's changes.
  * Changed the paths to local ../interface rather than global Higgs/ntuplepruducer/interface

2011 Sep 30
  Tagged with V00-05 - Wjets and Data with proper json files are produced with it.
  Attempted to properly include the TriggerObject - failed. Need to figure out the way to call the method from HLT trigger.
  Tested skimZ - works fine. 

2011 Sep 21
  Added the MET significance method to TCMET class. 

2011 aug 30
  Added updated json files for future productions
  Updated crab config files (still needs to be checked though)
  
  

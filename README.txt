Several packages need to be check out for the ntuple code to run properly.
This is true as of CMSSW_4_2_4.  

--- Ecal noise cleaning

cvs co -r V110523_BE -d DataFormats/AnomalousEcalDataFormats UserCode/csander/DataFormats/AnomalousEcalDataFormats
cvs co -r V110523_BE -d PhysicsTools/EcalAnomalousEventFilter UserCode/csander/PhysicsTools/EcalAnomalousEventFilter
cvs co -r V110523_BE -d Sandbox UserCode/csander/Sandbox
cvs co -r V19MAY2011_v3 JetMETAnalysis/ecalDeadCellTools


--- MET corrections

cvs co -r V04-04-04 JetMETCorrections/Type1MET


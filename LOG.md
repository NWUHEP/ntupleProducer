 Log of changes 
---------------

2013-Sep-17
 * Updated the code to work with CMSSW_5_3_11_patch6
 -Andrey

2013-Aug-18
 * Added back a few of the minor changes from Aug-01.  It was a waste of my time.  skimLepton was not even fully implemented into the analyzer,
 it was just floating around unused and uninitialized.  VarParser is still removed, looks like its more complicated then useful in
 its v6.4 iteration.  Trigger names still need to be readded for H->llgamma
 

2013-Aug-01
 * Using VarParser for selection isRealData option. Now one can run it with: ```cmsRun isRealData=1```
without editing the config file. One can also specify other option, or add (register) more.
 * Added an option for skimLepton from a cfg
   * Lowered the cut on electron pt skimming from 10 to 5
 * More HLT trigger names added for H -> llgamma dalits analysis.

2013-Jul-13
 * Ported this thing to git.

2012-Oct :  2013-Jul
 * Nobody likes to write the logs...

2012-Oct-04
 * Branch for 533 release

2011-Sep-30

 * Tagged with V00-05 - Wjets and Data with proper json files are produced with it.
 * Attempted to properly include the TriggerObject - failed. Need to figure out the way to call the method from HLT trigger.
 * Tested skimZ - works fine. 

2011-Sep-21
 * Added the MET significance method to TCMET class. 

2011-Aug-30
 * Added updated json files for future productions
 * Updated crab config files (still needs to be checked though)
  
  

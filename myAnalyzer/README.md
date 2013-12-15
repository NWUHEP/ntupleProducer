How to write an analyzer based on Northwestern ntuples
-----

Now have you produced the ntuples with our state-of-the-art ntuple Producer, what do you do next? How to analyze those?
Don't worry, we have you covered and this page will explain it all!

Presumably you got an output **nuTuple.root** file on hands.

* First, let's open it in Root:
```
[~] root -l nuTuple.root
```
You will get various _'no-dictionary'_ warnings. That's Okay, ignore them for now.


 * Next, in the root file **cd()** to a subdirectory where our main tree is saved and call **MakeSelector()** from that tree (you can browse the root file with TBrowser to find out that the tree we are interesed in is located under **ntupleProducer/eventTree**):
```
root [1] _file0.cd("ntupleProducer")
root [2] eventTree.MakeSelector("myAnalyzer")
```
That last command should give you something like:
```
Info in <TTreePlayer::MakeClass>: Files: myAnalyzer.h and myAnalyzer.C generated from TTree: eventTree
```
As you may have guessed, we will use TSelector method to write the analysis. More details on which can be found here: [http://root.cern.ch/root/html/TSelector.html][TSelector]
One may use other analysis methods, that's fine too, remember we have a rather simple root tree to deal with.

Two files have been created for us: **myAnalyzer.cc** and **myAnalyzer.h**

 * Now we need to edit **.h** file to let it know about TClone classes that we have inside our root file, just include those:
```cpp
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
```

 * For **.cc** file add the following into your **Process()** method:
```cpp
Bool_t myAnalyzer::Process(Long64_t entry)
{
  GetEntry(entry);
  for (Int_t i = 0; i < recoPhotons->GetSize(); ++i) {
    TCPhoton* thisPhoton = (TCPhoton*) recoPhotons->At(i);
    std::cout<<"Event number="<<eventNumber<<";  photon #"<<i<<", pt="<<thisPhoton->Pt()<<std::endl;
  }
```
That will get an entry and then loop over the collection of TCPhotons (it is called recoProtons in our tree).

 * How to run?
As it says in the Root's manual for TSelector, you run it with **tree->Process("myAnalyzer.c+")**.
However, that would work for a simple tree while ours contains various classes that root doesn't know about.
We need to load those before running and then do the Process().
Have a look at **run.C** script which contains all that is needed:
```cpp
void run(){
  gROOT->SetMacroPath(".:../src/:../interface/");

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

  fChain->Process("myAnalyzer.C+", "options");
}
```

 * You can run it simply as:  
```root -l -q run.C```
There is a python version as well: ```python run.py```

That's it!
Now you are ready to write your own analysis code. Here are other tips and tricks you may need to know, see below

Tips and tricks
---------------

#### How to fill the histograms

 * Define it in the begining of the myAnalyzer.cc
```cpp
TFile *f;
TH1F *h;
```
 * Initialize inside  myAnalyzer::Begin method:
```cpp
   f = new TFile("output.root","recreate");
   h = new TH1F("photon_pt","photon_pt", 100,0,100);
```

 * Fill it, somewhere in the Process():
```
h->Fill(thisPhoton->Pt());
```

 * Write and close the file in when Terminate()
```
f->Write();
f->Close();
```

In case one has to make a lot of histograms, one could use famous HistoManager,



#### How do I know which objects are available in the tree

 * Look into .h files for the class definitions.
 * Browse the root file with TBrowser, to see what's in there.

#### How do I know wjich methods are available fro each object:

You have to go and look into a paticulr classs where the object is defeined.
For example all methods available from Muons can be viewd
at ntupleProducer/sr/TCMuon.cc

#### How do I know that the method I call returns the value I want?
In order to make sure that you get what you want you need to go and find out how the particular object was filled and the varibles Set.
For example, in Muons we set the number of pixel hits
with _SetNumberOfValidPixelHits()_ method.
If you need to know where it's set you go to the **ntupleProducer.cc** code and find:
[Line 416][code-L416]

Most of the obcects we save in the ntuples are inherited from **TCPhysObject** class which itself is inherited from ROOT's TObject.
**TCPhoton** and **TCElectron** are inherited from a common **TCEgamma** class.

#### How to pass options to your analyzer

#### How to get the total number of events

#### How to use the trigger names

#### How to make easy histogramms with HistManager

#### How to use plugins

#### Do I need CMSSW to run it
 No, you only need Root.

[TSelector]: http://root.cern.ch/root/html/TSelector.html
[code-L416]: https://github.com/NWUHEP/ntupleProducer/blob/master/src/ntupleProducer.cc#L416

#include "Higgs/ntupleProducer/interface/ntupleProducer.h"

ntupleProducer::ntupleProducer(const edm::ParameterSet& iConfig)
{
    jetTag_           = iConfig.getUntrackedParameter<edm::InputTag>("JetTag");
    jptTag_           = iConfig.getUntrackedParameter<edm::InputTag>("JPTTag");
    metTag_           = iConfig.getUntrackedParameter<edm::InputTag>("METTag");
    metNoPUTag_       = iConfig.getUntrackedParameter<edm::InputTag>("METNoPUTag");
    muonTag_          = iConfig.getUntrackedParameter<edm::InputTag>("MuonTag");
    electronTag_      = iConfig.getUntrackedParameter<edm::InputTag>("ElectronTag");
    photonTag_        = iConfig.getUntrackedParameter<edm::InputTag>("PhotonTag");
    tauTag_           = iConfig.getUntrackedParameter<edm::InputTag>("TauTag");
    genJetTag_        = iConfig.getUntrackedParameter<edm::InputTag>("GenJetTag");
    primaryVtxTag_    = iConfig.getUntrackedParameter<edm::InputTag>("PrimaryVtxTag");
    rhoCorrTag_       = iConfig.getUntrackedParameter<edm::InputTag>("rhoCorrTag");
    hlTriggerResults_ = iConfig.getUntrackedParameter<string>("HLTriggerResults","TriggerResults");
    hltProcess_       = iConfig.getUntrackedParameter<string>("hltName");
    triggerPaths_     = iConfig.getUntrackedParameter<vector<string> >("triggers");

    saveJets_         = iConfig.getUntrackedParameter<bool>("saveJets");
    saveElectrons_    = iConfig.getUntrackedParameter<bool>("saveElectrons");
    saveMuons_        = iConfig.getUntrackedParameter<bool>("saveMuons");
    saveTaus_         = iConfig.getUntrackedParameter<bool>("saveTaus");
    savePhotons_      = iConfig.getUntrackedParameter<bool>("savePhotons");
    saveMET_          = iConfig.getUntrackedParameter<bool>("saveMET");
    saveGenJets_      = iConfig.getUntrackedParameter<bool>("saveGenJets");

    ecalFilterTag_    = iConfig.getUntrackedParameter<edm::InputTag>("ecalFilterTag");
    hcalFilterTag_    = iConfig.getUntrackedParameter<edm::InputTag>("hcalFilterTag");

    photonIsoCalcTag_ = iConfig.getParameter<edm::ParameterSet>("photonIsoCalcTag");
}

ntupleProducer::~ntupleProducer()
{

}

//
// member functions
//

// ------------ method called to for each event  ------------
void ntupleProducer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    eventNumber  = iEvent.id().event(); 
    runNumber    = iEvent.id().run();
    lumiSection  = (unsigned int)iEvent.getLuminosityBlock().luminosityBlock();
    bunchCross   = (unsigned int)iEvent.bunchCrossing();
    isRealData   = iEvent.isRealData();

    edm::Handle<reco::BeamSpot> beamSpotHandle;
    iEvent.getByLabel("offlineBeamSpot", beamSpotHandle);
    reco::BeamSpot vertexBeamSpot = *beamSpotHandle;

    beamSpot->SetXYZ(vertexBeamSpot.x0(), vertexBeamSpot.y0(), vertexBeamSpot.z0());

    int vtxCount, jetCount, jptCount, metCount, muCount, eleCount, photonCount, tauCount, genCount, genPartCount;
    vtxCount = jetCount = jptCount = metCount = muCount = eleCount = photonCount = tauCount = genCount = genPartCount = 0;
    float primaryVertexZ = -999;


    //////////////////////////
    //Get vertex information//
    //////////////////////////

    Handle<reco::VertexCollection> primaryVtcs;
    iEvent.getByLabel(primaryVtxTag_, primaryVtcs);

    for(VertexCollection::const_iterator iVtx = primaryVtcs->begin(); iVtx!= primaryVtcs->end(); ++iVtx){
        reco::Vertex myVtx = reco::Vertex(*iVtx);
        if(!myVtx.isValid() || myVtx.isFake()) continue;
        TCPrimaryVtx* vtxCon = new ((*primaryVtx)[vtxCount]) TCPrimaryVtx;
        vtxCon->SetPosition(myVtx.x(), myVtx.y(), myVtx.z());
        vtxCon->SetNDof(myVtx.ndof());
        vtxCon->SetChi2(myVtx.chi2());
        vtxCon->SetNtracks(myVtx.nTracks()); 
        vtxCon->SetSumPt2Trks(sumPtSquared(myVtx));
        if(vtxCount == 0) primaryVertexZ = myVtx.z();
        ++vtxCount;
    }


    ///////////////////////
    //get jet information//
    ///////////////////////

    Handle<double> rhoCorr;
    iEvent.getByLabel(rhoCorrTag_, rhoCorr);
    rhoFactor = (float)(*rhoCorr);

    if(saveJets_){

        Handle<vector<pat::Jet> > jets;
        iEvent.getByLabel(jetTag_, jets);

        for (vector<pat::Jet>::const_iterator iJet = jets->begin(); iJet != jets->end(); ++iJet) {

            if (iJet->pt() < 10.) continue;

            TCJet* jetCon = new ((*recoJets)[jetCount]) TCJet;

            //cout << "Uncorrected jet pt: " << iJet->correctedJet(0).pt() 
            //    << ", corrected jet pt: " << iJet->correctedJet(3).pt() 
            //    << endl; 

            jetCon->SetP4(iJet->correctedJet(0).px(), iJet->correctedJet(0).py(), iJet->correctedJet(0).pz(), iJet->correctedJet(0).energy());
            jetCon->SetVtx(0., 0., 0.);

            //cout << iJet->getSpecific.Zch << endl;;

            jetCon->SetChHadFrac(iJet->chargedHadronEnergyFraction());
            jetCon->SetNeuHadFrac(iJet->neutralHadronEnergyFraction());
            jetCon->SetChEmFrac(iJet->chargedEmEnergyFraction());
            jetCon->SetNeuEmFrac(iJet->neutralEmEnergyFraction());
            jetCon->SetNumConstit(iJet->chargedMultiplicity() + iJet->neutralMultiplicity());
            jetCon->SetNumChPart(iJet->chargedMultiplicity());

            jetCon->SetBDiscrTCHE(iJet->bDiscriminator("trackCountingHighEffBJetTags"));
            jetCon->SetBDiscrTCHP(iJet->bDiscriminator("trackCountingHighPurBJetTags"));
            jetCon->SetBDiscrSSVHE(iJet->bDiscriminator("simpleSecondaryVertexHighEffBJetTags"));
            jetCon->SetBDiscrJBP(iJet->bDiscriminator("jetProbabilityBJetTags"));

            jetCon->SetJetCorr(1, iJet->jecFactor(1));
            jetCon->SetJetCorr(2, iJet->jecFactor(2));
            jetCon->SetJetCorr(3, iJet->jecFactor(3));


            if (isRealData) {
                jetCon->SetJetCorr(4, iJet->jecFactor(4));
            } else {
                jetCon->SetJetCorr(4, 1.);
                jetCon->SetJetFlavor(iJet->partonFlavour());
            }

            /////////////////////////
            // Associate to vertex //
            /////////////////////////

            if(fabs(iJet->eta()) < 2.5){
                associateJetToVertex(*iJet, primaryVtcs, jetCon);
            } else {
                jetCon->SetVtxSumPtFrac(-1);
                jetCon->SetVtxSumPt(-1);
                jetCon->SetVtxTrackFrac(-1);
                jetCon->SetVtxNTracks(-1);
                jetCon->SetVtxSumPtIndex(0);
                jetCon->SetVtxCountIndex(0);
            }
            ++jetCount;
        }   
    }

    /////////////
    // Get MET //
    /////////////

    if (saveMET_){


        Handle<vector<pat::MET> > MET;
        iEvent.getByLabel(metTag_, MET);
        vector<pat::MET>::const_iterator met = MET->begin();

        if (MET->begin() != MET->end()) {
            recoMET->SetSumEt(met->sumEt());
            recoMET->SetMet(met->et());
            recoMET->SetPhi(met->phi());
            recoMET->SetMuonFraction(met->MuonEtFraction());
            recoMET->SetNeutralHadronFraction(met->NeutralHadEtFraction());
            recoMET->SetNeutralEMFraction(met->NeutralEMFraction());
            recoMET->SetChargedHadronFraction(met->ChargedHadEtFraction());
            recoMET->SetChargedEMFraction(met->ChargedEMEtFraction());

            //cout << met->nCorrections() << endl;

            //recoMET->SetCorrectedSumEt(iMET.sumEt());
            //recoMET->SetCorrectedMet(iMET.et());
            //recoMET->SetCorrectedPhi(iMET.phi());
        }

        Handle<vector<pat::MET> > METNoPU;
        iEvent.getByLabel(metNoPUTag_, METNoPU);
        met = METNoPU->begin();

        if (METNoPU->begin() != METNoPU->end()) {
            recoMETNoPU->SetSumEt(met->sumEt());
            recoMETNoPU->SetMet(met->et());
            recoMETNoPU->SetPhi(met->phi());
            recoMETNoPU->SetMuonFraction(met->MuonEtFraction());
            recoMETNoPU->SetNeutralHadronFraction(met->NeutralHadEtFraction());
            recoMETNoPU->SetNeutralEMFraction(met->NeutralEMFraction());
            recoMETNoPU->SetChargedHadronFraction(met->ChargedHadEtFraction());
            recoMETNoPU->SetChargedEMFraction(met->ChargedEMEtFraction());

            //cout << met->nCorrections() << endl;

            //recoMETNoPU->SetCorrectedSumEt(iMET.sumEt());
            //recoMETNoPU->SetCorrectedMet(iMET.et());
            //recoMETNoPU->SetCorrectedPhi(iMET.phi());
        }
    }

    ///////////////
    // Get muons //
    ///////////////

    if (saveMuons_) {

        //Handle<vector<pat::Muon> > muons;
        Handle<vector<reco::Muon> > muons;
        iEvent.getByLabel(muonTag_, muons);

        //for (vector<pat::Muon>::const_iterator iMuon = muons->begin(); iMuon != muons->end(); ++iMuon) {
        for (vector<reco::Muon>::const_iterator iMuon = muons->begin(); iMuon != muons->end(); ++iMuon) {
            if (!(iMuon->isGlobalMuon() && iMuon->isTrackerMuon()) || iMuon->pt() < 10.) continue;

            TCMuon* muCon = new ((*recoMuons)[muCount]) TCMuon;

            muCon->SetP4(iMuon->px(), iMuon->py(), iMuon->pz(), iMuon->energy());
            muCon->SetVtx(iMuon->globalTrack()->vx(),iMuon->globalTrack()->vy(),iMuon->globalTrack()->vz());
            muCon->SetPtError(iMuon->globalTrack()->ptError());
            muCon->SetCharge(iMuon->charge());
            muCon->SetIsGLB(iMuon->isGlobalMuon());
            muCon->SetIsTRK(iMuon->isTrackerMuon());
            muCon->SetNumberOfMatches(iMuon->numberOfMatches());
            muCon->SetNumberOfValidPixelHits(iMuon->globalTrack()->hitPattern().numberOfValidPixelHits());
            muCon->SetNumberOfValidTrackerHits(iMuon->globalTrack()->hitPattern().numberOfValidTrackerHits()); 
            muCon->SetNumberOfValidMuonHits(iMuon->globalTrack()->hitPattern().numberOfValidMuonHits());
            muCon->SetNumberOfLostPixelHits(iMuon->globalTrack()->hitPattern().numberOfLostPixelHits());
            muCon->SetNumberOfLostTrackerHits(iMuon->globalTrack()->hitPattern().numberOfLostTrackerHits());
            muCon->SetNormalizedChi2(iMuon->globalTrack()->normalizedChi2());

            muCon->SetCaloComp(iMuon->caloCompatibility());
            muCon->SetSegComp(muon::segmentCompatibility(*iMuon));

            muCon->SetNtracks03(iMuon->isolationR03().nTracks);
            muCon->SetEmIso03(iMuon->isolationR03().emEt);
            muCon->SetHadIso03(iMuon->isolationR03().hadEt);
            muCon->SetTrkIso03(iMuon->isolationR03().sumPt);

            muCon->SetNtracks05(iMuon->isolationR05().nTracks);
            muCon->SetEmIso05(iMuon->isolationR05().emEt);
            muCon->SetHadIso05(iMuon->isolationR05().hadEt);
            muCon->SetTrkIso05(iMuon->isolationR05().sumPt);

            muCon->SetPfSumPt(0.3, iMuon->pfIsolationR03().sumChargedParticlePt);
            muCon->SetPfEGamma(0.3, iMuon->pfIsolationR03().sumPhotonEt);
            muCon->SetPfENeutral(0.3, iMuon->pfIsolationR03().sumNeutralHadronEt);
            muCon->SetPfSumPt(0.4, iMuon->pfIsolationR04().sumChargedParticlePt);
            muCon->SetPfEGamma(0.4, iMuon->pfIsolationR04().sumPhotonEt);
            muCon->SetPfENeutral(0.4, iMuon->pfIsolationR04().sumNeutralHadronEt);

            muCount++;
        }
    }


    ///////////////////
    // Get electrons //
    ///////////////////

    if (saveElectrons_) {

        Handle<vector<pat::Electron> > electrons;
        iEvent.getByLabel(electronTag_, electrons);

        for (vector<pat::Electron>::const_iterator iElectron = electrons->begin(); iElectron != electrons->end(); ++iElectron) {
            if (iElectron->pt() < 10) continue;

            TCElectron* eleCon = new ((*recoElectrons)[eleCount]) TCElectron;

            eleCon->SetP4(iElectron->px(),iElectron->py(),iElectron->pz(),iElectron->p());
            eleCon->SetVtx(iElectron->gsfTrack()->vx(),iElectron->gsfTrack()->vy(),iElectron->gsfTrack()->vz());
            eleCon->SetCharge(iElectron->charge());

            eleCon->SetNumberOfValidPixelHits(iElectron->gsfTrack()->hitPattern().numberOfValidPixelHits());
            eleCon->SetNumberOfValidTrackerHits(iElectron->gsfTrack()->hitPattern().numberOfValidTrackerHits());
            eleCon->SetNumberOfLostPixelHits(iElectron->gsfTrack()->hitPattern().numberOfLostPixelHits());
            eleCon->SetNumberOfLostTrackerHits(iElectron->gsfTrack()->hitPattern().numberOfLostTrackerHits());

            eleCon->SetIsEB(iElectron->isEB());
            eleCon->SetIsEE(iElectron->isEE());
            eleCon->SetIsInGap(iElectron->isGap());

            eleCon->SetEmIso03(iElectron->dr03EcalRecHitSumEt());
            eleCon->SetHadIso03(iElectron->dr03HcalTowerSumEt());
            eleCon->SetTrkIso03(iElectron->dr03TkSumPt());
            eleCon->SetEmIso04(iElectron->dr04EcalRecHitSumEt());
            eleCon->SetHadIso04(iElectron->dr04HcalTowerSumEt());
            eleCon->SetTrkIso04(iElectron->dr04TkSumPt());

            eleCon->SetHadOverEm(iElectron->hadronicOverEm());
            eleCon->SetDphiSuperCluster(iElectron->deltaPhiSuperClusterTrackAtVtx());
            eleCon->SetDetaSuperCluster(iElectron->deltaEtaSuperClusterTrackAtVtx());
            eleCon->SetSigmaIetaIeta(iElectron->sigmaIetaIeta());
            eleCon->SetFBrem(iElectron->fbrem());
            eleCon->SetEOverP(iElectron->eSuperClusterOverP());

            eleCon->SetConversionFlag(iElectron->convFlags());
            eleCon->SetConversionDist(iElectron->convDist());
            eleCon->SetConversionDcot(iElectron->convDcot());
            eleCon->SetConversionRad(iElectron->convRadius());

            eleCon->SetCutLevel(iElectron->electronID("eidVBTF95"), 95);
            eleCon->SetCutLevel(iElectron->electronID("eidVBTF90"), 90);
            eleCon->SetCutLevel(iElectron->electronID("eidVBTF85"), 85);
            eleCon->SetCutLevel(iElectron->electronID("eidVBTF80"), 80);
            eleCon->SetCutLevel(iElectron->electronID("eidVBTF70"), 70);
            eleCon->SetCutLevel(iElectron->electronID("eidVBTF60"), 60);

            eleCon->SetPfEGamma(0.3, iElectron->pfIsolationVariables().photonIso);
            eleCon->SetPfSumPt(0.3, iElectron->pfIsolationVariables().chargedHadronIso);
            eleCon->SetPfENeutral(0.3, iElectron->pfIsolationVariables().neutralHadronIso);

            eleCount++;
        }
    }

    /////////////////
    // Get photons //
    /////////////////

    if (savePhotons_) {
        Handle<vector<pat::Photon> > photons;
        iEvent.getByLabel(photonTag_, photons);

        for (vector<pat::Photon>::const_iterator iPhoton = photons->begin(); iPhoton != photons->end() ; ++iPhoton) {

            TCPhoton* myPhoton = new ((*recoPhotons)[photonCount]) TCPhoton();
            myPhoton->SetP4(iPhoton->px(), iPhoton->py(), iPhoton->pz(), iPhoton->p());
            myPhoton->SetVtx(iPhoton->vx(), iPhoton->vy(), iPhoton->vz());
            myPhoton->SetHadOverEm(iPhoton->hadronicOverEm());
            myPhoton->SetSigmaIEtaIEta(iPhoton->sigmaIetaIeta());
            myPhoton->SetR9(iPhoton->r9());
            myPhoton->SetEtaSupercluster(iPhoton->superCluster()->eta());
            myPhoton->SetTrackVeto(iPhoton->hasPixelSeed());

            //Conversion info
            reco::ConversionRefVector conversions = iPhoton->conversions();
            int   conversionCount = 0;

            for (reco::ConversionRefVector::const_iterator iConversion = conversions.begin(); iConversion != conversions.end(); ++iConversion) {
                const reco::ConversionRef myConversion = *iConversion;
                ++conversionCount;
            }
            myPhoton->SetNumberOfConversions(conversionCount);

            PhotonIsolationCalculator* photonIsolationCalculator = new PhotonIsolationCalculator();
            photonIsolationCalculator->setup(photonIsoCalcTag_);
            reco::Photon::FiducialFlags fiducialFlags;
            reco::Photon::IsolationVariables isolVarR03, isolVarR04;
            int vtxPhoCounter = 0;

            for(VertexCollection::const_iterator iVtx = primaryVtcs->begin(); iVtx!= primaryVtcs->end(); ++iVtx){
                float *isoVtxAB;
                reco::Vertex myVtx = reco::Vertex(*iVtx);
                if(!myVtx.isValid() || myVtx.isFake()) continue;
                // Calculate fiducial flags and isolation variable. Blocked are filled from the isolationCalculator
                isoVtxAB = photonIsolationCalculator->calculateVtx( &(*iPhoton),iEvent,iSetup,fiducialFlags,isolVarR03, isolVarR04, myVtx);
                myPhoton->SetTRKIsoVtxDR03(isoVtxAB[1]);
                myPhoton->SetTRKIsoVtxDR04(isoVtxAB[0]);
                vtxPhoCounter++;

            }
            ++photonCount;
        }
    }

    //////////////
    // Get taus //
    //////////////

    if (saveTaus_) {
        Handle<vector<pat::Tau> > taus;
        iEvent.getByLabel(tauTag_, taus);

        for (vector<pat::Tau>::const_iterator iTau = taus->begin(); iTau != taus->end(); ++iTau) {            

            if (! (iTau->isPFTau()) ) continue;


            // reject taus with no charged track at the beginning of the loop
            if (iTau->signalPFChargedHadrCands().size() < 1) continue;

            if (iTau->pt() < 10) continue;

            TCTau* tauCon = new ((*recoTaus)[tauCount]) TCTau;

            tauCon->SetNChHad(iTau->signalPFChargedHadrCands().size());
            tauCon->SetNGamma (iTau->signalPFGammaCands().size());
            tauCon->SetNNeutrHad (iTau->signalPFNeutrHadrCands().size());
            tauCon->SetCharge(iTau->charge());
            tauCon->SetDecayMode(iTau->decayMode());

            tauCon->SetP4(iTau->px(),iTau->py(),iTau->pz(),iTau->energy());
            tauCon->SetCharge(iTau->charge());

            if (iTau->leadPFChargedHadrCand()->trackRef().isNonnull()) {
                tauCon->SetLeadChHadP4(iTau->leadPFChargedHadrCand()->px(),
                        iTau->leadPFChargedHadrCand()->py(),
                        iTau->leadPFChargedHadrCand()->pz(),
                        iTau->leadPFChargedHadrCand()->energy());

                tauCon->SetPositionFromTrack(iTau->leadPFChargedHadrCand()->trackRef()->vx(),
                        iTau->leadPFChargedHadrCand()->trackRef()->vy(),
                        iTau->leadPFChargedHadrCand()->trackRef()->vz());
            }


            if (iTau->signalPFGammaCands().size()+iTau->signalPFNeutrHadrCands().size()>0) 
                tauCon->SetLeadNeutrP4(iTau->leadPFNeutralCand()->px(),
                        iTau->leadPFNeutralCand()->py(),
                        iTau->leadPFNeutralCand()->pz(),
                        iTau->leadPFNeutralCand()->energy());

            tauCon->SetPositionFromTau(iTau->vx(),iTau->vy(), iTau->vz());

            tauCon->SetIsoGammaEtSum(iTau->isolationPFGammaCandsEtSum());
            tauCon->SetIsoChHadPtSum(iTau->isolationPFChargedHadrCandsPtSum());



            // set the discriminators
            // note that the strings for PAT and RECO are different. The names of the TCTau accessors are set following the RECO names
            // the "mapping" is taken from  tauTools.py 

            tauCon->SetHpsPFTauDiscriminationByDecayModeFinding(iTau->tauID("decayModeFinding"));  // "DiscriminationByDecayModeFinding"

            // isolation
            //
            tauCon->SetHpsPFTauDiscriminationByVLooseIsolation(iTau->tauID("byVLooseIsolation")); // "DiscriminationByVLooseIsolation"
            tauCon->SetHpsPFTauDiscriminationByLooseIsolation (iTau->tauID("byLooseIsolation"));  // "DiscriminationByLooseIsolation"
            tauCon->SetHpsPFTauDiscriminationByMediumIsolation(iTau->tauID("byMediumIsolation")); // "DiscriminationByMediumIsolation"  
            tauCon->SetHpsPFTauDiscriminationByTightIsolation (iTau->tauID("byTightIsolation"));  // "DiscriminationByTightIsolation"
            // isolation with corrections
            tauCon->SetHpsPFTauDiscriminationByVLooseIsolationDBSumPtCorr	     
                (iTau->tauID("byVLooseIsolationDeltaBetaCorr")); // "DiscriminationByVLooseIsolationDBSumPtCorr"
            tauCon->SetHpsPFTauDiscriminationByLooseIsolationDBSumPtCorr	     
                (iTau->tauID("byLooseIsolationDeltaBetaCorr")); // "DiscriminationByLooseIsolationDBSumPtCorr"
            tauCon->SetHpsPFTauDiscriminationByMediumIsolationDBSumPtCorr
                (iTau->tauID("byMediumIsolationDeltaBetaCorr")); // "DiscriminationByMediumIsolationDBSumPtCorr"
            tauCon->SetHpsPFTauDiscriminationByTightIsolationDBSumPtCorr
                (iTau->tauID("byTightIsolationDeltaBetaCorr")); // "DiscriminationByTightIsolationDBSumPtCorr"
            // combined isolation with corrections
            tauCon->SetHpsPFTauDiscriminationByVLooseCombinedIsolationDBSumPtCorr
                (iTau->tauID("byVLooseCombinedIsolationDeltaBetaCorr")); // "DiscriminationByVLooseCombinedIsolationDBSumPtCorr"
            tauCon->SetHpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr
                (iTau->tauID("byLooseCombinedIsolationDeltaBetaCorr")); // "DiscriminationByLooseCombinedIsolationDBSumPtCorr"
            tauCon->SetHpsPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr
                (iTau->tauID("byMediumCombinedIsolationDeltaBetaCorr")); // "DiscriminationByMediumCombinedIsolationDBSumPtCorr"
            tauCon->SetHpsPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr 
                (iTau->tauID("byTightCombinedIsolationDeltaBetaCorr")); // "DiscriminationByTightCombinedIsolationDBSumPtCorr"

            // anti e/mu discriminators
            tauCon->SetHpsPFTauDiscriminationAgainstElectronLoose (iTau->tauID("againstElectronLoose")); // "DiscriminationByLooseElectronRejection"
            tauCon->SetHpsPFTauDiscriminationAgainstElectronMedium(iTau->tauID("againstElectronMedium")); // "DiscriminationByMediumElectronRejection"
            tauCon->SetHpsPFTauDiscriminationAgainstElectronTight (iTau->tauID("againstElectronTight")); // "DiscriminationByTightElectronRejection"

            tauCon->SetHpsPFTauDiscriminationAgainstMuonLoose  (iTau->tauID("againstMuonLoose")); // "DiscriminationByLooseMuonRejection")
            //	  tauCon->SetHpsPFTauDiscriminationAgainstMuonMediumt(iTau->tauID("againstMuonMedium")); // "DiscriminationByMediumMuonRejection" <- not in python
            tauCon->SetHpsPFTauDiscriminationAgainstMuonTight  (iTau->tauID("againstMuonTight")); // "DiscriminationByTightMuonRejection"

            tauCount++;
        } // loop over taus
    } // if save taus



    ////////////////////////
    // Get gen-level info //
    ////////////////////////


    if (!isRealData) {

        Handle<GenEventInfoProduct> GenEventInfoHandle;
        iEvent.getByLabel("generator", GenEventInfoHandle);

        evtWeight = ptHat = qScale = -1;

        if (GenEventInfoHandle.isValid()) {
            //qScale       = GenEventInfoHandle->qScale();
            ptHat        = (GenEventInfoHandle->hasBinningValues() ? GenEventInfoHandle->binningValues()[0] : 0.0);
            //evtWeight    = GenEventInfoHandle->weight();
        }


        ////////////////////
        // PU information //
        ////////////////////


        Handle<std::vector< PileupSummaryInfo > >  PUInfo;
        iEvent.getByLabel(edm::InputTag("addPileupInfo"), PUInfo);
        std::vector<PileupSummaryInfo>::const_iterator iPV;

        for(iPV = PUInfo->begin(); iPV != PUInfo->end(); ++iPV){
          if (iPV->getBunchCrossing() == 0){
            nPUVertices     = iPV->getPU_NumInteractions();
            nPUVerticesTrue = iPV->getTrueNumInteractions();
          }
        }


        //////////////////////
        // Get genParticles //
        //////////////////////

        Handle<GenParticleCollection> genParticleColl;
        iEvent.getByLabel("genParticles", genParticleColl);

        for (GenParticleCollection::const_iterator iGenPart = genParticleColl->begin(); iGenPart != genParticleColl->end(); ++iGenPart) {
            const reco::GenParticle myParticle = reco::GenParticle(*iGenPart);

            if (myParticle.status() == 1 && ((abs(myParticle.pdgId()) >= 11 && abs(myParticle.pdgId()) <= 16) || myParticle.pdgId() == 22)) {
                TCGenParticle* genCon = new ((*genParticles)[genPartCount]) TCGenParticle;

                genCon->SetPosition(myParticle.vx(), myParticle.vy(), myParticle.vz() );
                genCon->SetP4(myParticle.px(), myParticle.py(), myParticle.pz(), myParticle.energy() );
                genCon->SetCharge(myParticle.charge());
                genCon->SetPDGId(myParticle.pdgId());
                genCon->SetMother(myParticle.mother()->pdgId());
                if (myParticle.pdgId() == 22 && myParticle.mother()->numberOfMothers() != 0) genCon->SetGrandmother(myParticle.mother()->mother()->pdgId());
                ++genPartCount;
            }


            if (myParticle.status() == 3 && (abs(myParticle.pdgId()) == 6 || abs(myParticle.pdgId()) == 39 || abs(myParticle.pdgId()) == 23 || abs(myParticle.pdgId()) == 24)) {
                for (size_t i = 0; i < myParticle.numberOfDaughters(); ++i) {
                    const reco::Candidate *myDaughter = myParticle.daughter(i);
                    if (myDaughter->status() != 1 && (abs(myDaughter->pdgId()) == 5 || (abs(myDaughter->pdgId()) >= 11 && abs(myDaughter->pdgId()) <= 16))) {
                        TCGenParticle* genCon = new ((*genParticles)[genPartCount]) TCGenParticle;

                        genCon->SetPosition(myDaughter->vx(), myDaughter->vy(), myDaughter->vz() );
                        genCon->SetP4(myDaughter->px(), myDaughter->py(), myDaughter->pz(), myDaughter->energy() );
                        genCon->SetCharge(myDaughter->charge());
                        genCon->SetPDGId(myDaughter->pdgId());
                        genCon->SetMother(myParticle.pdgId());
                        ++genPartCount;
                    }
                }
            }
        }


        /////////////////
        // Get genJets //
        /////////////////

        if (saveGenJets_) {

            Handle<reco::GenJetCollection> GenJets;
            iEvent.getByLabel(genJetTag_, GenJets);

            for (GenJetCollection::const_iterator iJet = GenJets->begin(); iJet!= GenJets->end(); ++iJet) {
                reco::GenJet myJet = reco::GenJet(*iJet);      
                if (myJet.pt() > 10) { 
                    TCGenJet* jetCon = new ((*genJets)[genCount]) TCGenJet;
                    jetCon->SetP4(myJet.px(), myJet.py(), myJet.pz(), myJet.energy());
                    jetCon->SetHadEnergy(myJet.hadEnergy());
                    jetCon->SetEmEnergy(myJet.emEnergy());
                    jetCon->SetInvEnergy(myJet.invisibleEnergy());
                    jetCon->SetAuxEnergy(myJet.auxiliaryEnergy());
                    jetCon->SetNumConstit(myJet.getGenConstituents().size());
                }
                ++genCount;	
            }
        }
    }


    ///////////////////
    // Noise filters //
    ///////////////////

    if (isRealData) {

        Handle<bool> hcalNoiseFilterHandle;
        iEvent.getByLabel(hcalFilterTag_, hcalNoiseFilterHandle);
        if (hcalNoiseFilterHandle.isValid())  isNoiseHcal = !(Bool_t)(*hcalNoiseFilterHandle);

        isDeadEcalCluster = kFALSE;
        Handle<AnomalousECALVariables> anomalousECALvarsHandle;
        iEvent.getByLabel(ecalFilterTag_, anomalousECALvarsHandle);
        AnomalousECALVariables anomalousECALvars;

        if (anomalousECALvarsHandle.isValid()) {
            anomalousECALvars = *anomalousECALvarsHandle;
            isDeadEcalCluster = anomalousECALvars.isDeadEcalCluster();
        }

        edm::Handle<BeamHaloSummary> TheBeamHaloSummary;
        iEvent.getByLabel("BeamHaloSummary",TheBeamHaloSummary);
        const BeamHaloSummary TheSummary = (*TheBeamHaloSummary.product() );

        isCSCTightHalo = TheSummary.CSCTightHaloId();
        isCSCLooseHalo = TheSummary.CSCLooseHaloId();

        isScraping = isFilteredOutScraping(iEvent, iSetup, 10, 0.25); 
    }

    ////////////////////////////  
    // get trigger information//
    ////////////////////////////

    edm::Handle<TriggerResults> hltR;
    triggerResultsTag_ = InputTag(hlTriggerResults_,"",hltProcess_);
    iEvent.getByLabel(triggerResultsTag_,hltR);

    const TriggerNames & triggerNames = iEvent.triggerNames(*hltR);
    hlNames=triggerNames.triggerNames();   

    triggerStatus = 0x0;    

    for (int i=0; i < (int)hlNames.size(); ++i) {      
        if (!triggerDecision(hltR, i)) continue;	
        for (int j = 0; j < (int)triggerPaths_.size(); ++j){
            if (triggerPaths_[j] == "") continue;
            if (hlNames[i].compare(0, triggerPaths_[j].length(),triggerPaths_[j]) == 0) {
                //cout << hlNames[i] << " ?= " << triggerPaths_[j] << endl;
                triggerStatus |= 0x01 << j;
                if (isRealData) {
                    pair<int, int> preScales;
                    preScales = hltConfig_.prescaleValues(iEvent, iSetup, hlNames[i]); 
                    hltPrescale[j] = preScales.first*preScales.second;
                }
            }
        }
    } 
    //printf("%016llX\n", triggerStatus);

    edm::Handle<trigger::TriggerEvent> triggerEvents;
    iEvent.getByLabel("hltTriggerSummaryAOD",triggerEvents);
    trigger::TriggerObjectCollection triggerObjCol = triggerEvents->getObjects();
    int triggerCount = 0;

    for(trigger::TriggerObjectCollection::const_iterator iTrigObj = triggerObjCol.begin(); iTrigObj != triggerObjCol.end(); ++iTrigObj) { 
        TCTriggerObject * thisTrig = new ((*triggerObjects)[triggerCount]) TCTriggerObject;
        thisTrig->setId(iTrigObj->id());
        thisTrig->setP4(iTrigObj->px(), iTrigObj->py(), iTrigObj->pz(), iTrigObj->energy());
        ++triggerCount;
    }

    ++nEvents;

    if (eleCount > 0 || muCount > 0) eventTree -> Fill(); // possibly specify a cut in configuration

    primaryVtx->Clear("C");
    recoJets->Clear("C");
    recoJPT->Clear("C");
    recoMuons->Clear("C");
    recoElectrons->Clear("C");
    recoTaus->Clear("C");
    recoPhotons->Clear("C");
    triggerObjects->Clear("C");
    genJets->Clear("C");
    genParticles->Clear("C");
}

// ------------ method called once each job just before starting event loop  ------------
void  ntupleProducer::beginJob()
{  
    eventTree      = fs->make<TTree>("eventTree","eventTree");
    runTree        = fs->make<TTree>("runTree","runTree");
    jobTree        = fs->make<TTree>("jobTree", "jobTree");

    primaryVtx     = new TClonesArray("TCPrimaryVtx");
    recoJets       = new TClonesArray("TCJet");
    recoJPT        = new TClonesArray("TCJet");
    recoElectrons  = new TClonesArray("TCElectron");
    recoMuons      = new TClonesArray("TCMuon");
    recoTaus       = new TClonesArray("TCTau");
    recoPhotons    = new TClonesArray("TCPhoton");
    triggerObjects = new TClonesArray("TCTriggerObject");
    genJets        = new TClonesArray("TCGenJet");
    genParticles   = new TClonesArray("TCGenParticle");
    beamSpot       = new TVector3();
    recoMET        = 0;
    recoMETNoPU    = 0;

    eventTree->Branch("recoJets",&recoJets, 6400, 0);
    eventTree->Branch("recoJPT",&recoJPT, 6400, 0);
    eventTree->Branch("recoElectrons",&recoElectrons, 6400, 0);
    eventTree->Branch("recoMuons",&recoMuons, 6400, 0);
    eventTree->Branch("recoTaus",&recoTaus, 6400, 0);
    eventTree->Branch("recoPhotons",&recoPhotons, 6400, 0);
    eventTree->Branch("recoMET", &recoMET, 6400, 0);
    eventTree->Branch("recoMETNoPU", &recoMETNoPU, 6400, 0);
    eventTree->Branch("triggerObjects", &triggerObjects, 6400, 0);
    eventTree->Branch("genJets",&genJets, 6400, 0);
    eventTree->Branch("genParticles",&genParticles, 6400, 0);

    eventTree->Branch("primaryVtx",&primaryVtx, 6400, 0);
    eventTree->Branch("beamSpot", &beamSpot, 6400, 0);
    eventTree->Branch("nPUVertices", &nPUVertices, "nPUVertices/I");
    eventTree->Branch("nPUVerticesTrue", &nPUVerticesTrue, "nPUVerticesTrue/F");

    eventTree->Branch("isRealData",&isRealData, "isRealData/O");
    eventTree->Branch("runNumber",&runNumber, "runNumber/l");
    eventTree->Branch("eventNumber",&eventNumber, "eventNumber/i");
    eventTree->Branch("lumiSection",&lumiSection, "lumiSection/e");
    eventTree->Branch("bunchCross",&bunchCross, "bunchCross/i");

    eventTree->Branch("isScraping",&isScraping, "isScraping/O");
    eventTree->Branch("isNoiseHcal",&isNoiseHcal, "isNoiseHcal/O");
    eventTree->Branch("isDeadEcalCluster",&isDeadEcalCluster, "isDeadEcalCluster/O");
    eventTree->Branch("isCSCTightHalo",&isCSCTightHalo, "isCSCTightHalo/O");
    eventTree->Branch("isCSCLooseHalo",&isCSCLooseHalo, "isCSCLooseHalo/O");

    eventTree->Branch("ptHat",&ptHat, "ptHat/F");
    eventTree->Branch("qScale", &qScale, "qScale/F");
    eventTree->Branch("evtWeight", &evtWeight, "evtWeight/F");
    eventTree->Branch("rhoFactor",&rhoFactor, "rhoFactor/F");
    eventTree->Branch("triggerStatus",&triggerStatus, "triggerStatus/l");
    eventTree->Branch("hltPrescale",hltPrescale, "hltPrescale[64]/i");

    runTree->Branch("deliveredLumi",&deliveredLumi, "deliveredLumi/F");
    runTree->Branch("recordedLumi",&recordedLumi, "recordedLumi/F");
    runTree->Branch("runNumber",&runNumber, "runNumber/l");

    jobTree->Branch("savedTriggerNames",savedTriggerNames, "savedTriggerNames[64]/C");
    jobTree->Branch("nEvents",&nEvents, "nEvents/i");

    // Initialize HLT prescales //

    for (int i = 0; i < (int)(sizeof(hltPrescale)/sizeof(int)); ++i) hltPrescale[i] = 1;

    // Start counting number of events per job //
    nEvents = 0;
}

void ntupleProducer::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup)
{
    bool changed = true; 
    hltConfig_.init(iRun, iSetup, hltProcess_, changed);
    deliveredLumi = 0;
    recordedLumi  = 0;
}

void ntupleProducer::endLuminosityBlock(const edm::LuminosityBlock& iLumi, const edm::EventSetup& iSetup)
{
    if (isRealData) {
        edm::Handle<LumiSummary> lumiSummary;
        iLumi.getByLabel("lumiProducer", lumiSummary);

        deliveredLumi  += lumiSummary->avgInsDelLumi()*93.244;
        recordedLumi   += deliveredLumi*lumiSummary->liveFrac();
    }
}

void ntupleProducer::endRun(const edm::Run& iRun, const edm::EventSetup& iSetup)
{
    //cout<<"\t Integrated luminosity = "<<deliveredLumi<<endl;
    runTree->Fill();
}

void ntupleProducer::endJob() 
{
    for (int i =0; i < (int)triggerPaths_.size(); ++i) savedTriggerNames[i] = triggerPaths_[i];
    cout<<nEvents<<endl;
    jobTree->Fill();
}

bool ntupleProducer::triggerDecision(edm::Handle<edm::TriggerResults> &hltR, int iTrigger)
{
    bool triggerPassed = false;
    if(hltR->wasrun(iTrigger) &&
            hltR->accept(iTrigger) &&
            !hltR->error(iTrigger) ){
        triggerPassed = true;
    }
    return triggerPassed;
}

float ntupleProducer::sumPtSquared(const Vertex & v)
{
    float sum = 0.;
    float pT;
    for (Vertex::trackRef_iterator it = v.tracks_begin(); it != v.tracks_end(); it++) {
        pT = (**it).pt();
        float epT=(**it).ptError(); pT=pT>epT ? pT-epT : 0;

        sum += pT*pT;
    }
    return sum;
}

void ntupleProducer::associateJetToVertex(pat::Jet inJet, Handle<reco::VertexCollection> vtxCollection, TCJet *outJet)
{
    const reco::TrackRefVector &tracks = inJet.associatedTracks(); 


    vector<float>  associatedTrackSumPt;
    vector<float>  associatedTrackCount;
    vector<const reco::Track*> jetTracks;
    float sumTrackX, sumTrackY, sumTrackZ, sumTrackPt;
    int   nJetTracks = 0;
    int   vCount = 0;

    sumTrackX = sumTrackY = sumTrackZ  = sumTrackPt = 0;

    for (TrackRefVector::const_iterator iTrack = tracks.begin(); iTrack != tracks.end(); ++iTrack) {
        const reco::Track &jetTrack = **iTrack;

        sumTrackPt += jetTrack.pt();
        sumTrackX  += jetTrack.vx();
        sumTrackY  += jetTrack.vy();            
        sumTrackZ  += jetTrack.vz();
        jetTracks.push_back(&jetTrack);
        ++nJetTracks;
    }

    if(jetTracks.size() == 0){
        outJet->SetVtxSumPtFrac(-1);
        outJet->SetVtxSumPt(0);
        outJet->SetVtxTrackFrac(-1);
        outJet->SetVtxNTracks(0);
        outJet->SetVtxSumPtIndex(0);
        outJet->SetVtxCountIndex(0);
        outJet->SetVtx(0., 0., 0.);      	
    } else {
        outJet->SetVtx(sumTrackX/nJetTracks, sumTrackY/nJetTracks, sumTrackZ/nJetTracks);       
        for (VertexCollection::const_iterator iVtx = vtxCollection->begin(); iVtx!= vtxCollection->end(); ++iVtx) {	      
            reco::Vertex myVtx = reco::Vertex(*iVtx); 
            if(!myVtx.isValid() || myVtx.isFake()) continue;
            associatedTrackSumPt.push_back(0);            
            associatedTrackCount.push_back(0);            

            for(Vertex::trackRef_iterator iTrackRef = myVtx.tracks_begin(); iTrackRef != myVtx.tracks_end(); ++iTrackRef){
                const edm::RefToBase<reco::Track> &myTrackRef = *iTrackRef; 

                if(myTrackRef.isAvailable()){
                    const reco::Track &myVertexTrack = *myTrackRef.get();		

                    for(vector<const reco::Track*>::const_iterator iTrack = jetTracks.begin(); iTrack != jetTracks.end(); ++iTrack){
                        if (*iTrack == &myVertexTrack) {
                            associatedTrackSumPt.at(vCount) += myVertexTrack.pt()/sumTrackPt; 
                            associatedTrackCount.at(vCount) += 1/nJetTracks; 
                        }
                    }
                }
            }
            ++vCount;  
        }

        float maxSumPtFraction = 0; float maxCountFraction = 0;
        int   vtxSumPtIndex = 0; int vtxCountIndex = 0;
        int count = 0;

        for (int i = 0; i < vCount; ++i) {
            if (associatedTrackSumPt.at(i) > maxSumPtFraction) {
                maxSumPtFraction = associatedTrackSumPt.at(i);   
                vtxSumPtIndex = count + 1;
            }
            if (associatedTrackCount.at(i) > maxCountFraction) {
                maxCountFraction = associatedTrackCount.at(i);   
                vtxCountIndex = count + 1;
            }
            ++count;
        }
        outJet->SetVtxSumPtFrac(maxSumPtFraction);
        outJet->SetVtxSumPt(sumTrackPt);
        outJet->SetVtxTrackFrac(maxCountFraction);
        outJet->SetVtxNTracks(nJetTracks);
        outJet->SetVtxSumPtIndex(vtxSumPtIndex);
        outJet->SetVtxCountIndex(vtxCountIndex);
    }

}


bool ntupleProducer::isFilteredOutScraping( const edm::Event& iEvent, const edm::EventSetup& iSetup, int numtrack, double thresh)
{

    bool  accepted = false;
    float fraction = 0;  
    // get GeneralTracks collection

    edm::Handle<reco::TrackCollection> tkRef;
    iEvent.getByLabel("generalTracks",tkRef);    
    const reco::TrackCollection* tkColl = tkRef.product();

    int numhighpurity=0;
    reco::TrackBase::TrackQuality _trackQuality = reco::TrackBase::qualityByName("highPurity");

    if(tkColl->size()>(UInt_t)numtrack){ 
        reco::TrackCollection::const_iterator itk = tkColl->begin();
        reco::TrackCollection::const_iterator itk_e = tkColl->end();
        for(;itk!=itk_e;++itk){
            if(itk->quality(_trackQuality)) numhighpurity++;
        }
        fraction = (float)numhighpurity/(float)tkColl->size();
        if(fraction>thresh) accepted=true;
    } else {
        //if less than 10 Tracks accept the event anyway    
        accepted= true;
    }
    return !accepted;  //if filtered out it's not accepted.
}

//define this as a plug-in
DEFINE_FWK_MODULE(ntupleProducer);

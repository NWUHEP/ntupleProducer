#include "../interface/ntupleProducer.h"
#include "DataFormats/Math/interface/deltaR.h"

ntupleProducer::ntupleProducer(const edm::ParameterSet& iConfig)
{
    jetTag_           = iConfig.getUntrackedParameter<edm::InputTag>("JetTag");
    metTag_           = iConfig.getUntrackedParameter<edm::InputTag>("METTag");
    muonTag_          = iConfig.getUntrackedParameter<edm::InputTag>("MuonTag");
    electronTag_      = iConfig.getUntrackedParameter<edm::InputTag>("ElectronTag");
    photonTag_        = iConfig.getUntrackedParameter<edm::InputTag>("PhotonTag");
    tauTag_           = iConfig.getUntrackedParameter<edm::InputTag>("TauTag");
    genJetTag_        = iConfig.getUntrackedParameter<edm::InputTag>("GenJetTag");
    primaryVtxTag_    = iConfig.getUntrackedParameter<edm::InputTag>("PrimaryVtxTag");

    rhoCorrTag_       = iConfig.getUntrackedParameter<edm::InputTag>("rhoCorrTag");
    rho25CorrTag_     = iConfig.getUntrackedParameter<edm::InputTag>("rho25CorrTag");
    hlTriggerResults_ = iConfig.getUntrackedParameter<string>("HLTriggerResults","TriggerResults");
    hltProcess_       = iConfig.getUntrackedParameter<string>("hltName");
    triggerPaths_     = iConfig.getUntrackedParameter<vector<string> >("triggers");

    partFlowTag_      = iConfig.getUntrackedParameter<edm::InputTag>("partFlowTag");

    saveJets_         = iConfig.getUntrackedParameter<bool>("saveJets");
    saveElectrons_    = iConfig.getUntrackedParameter<bool>("saveElectrons");
    saveMuons_        = iConfig.getUntrackedParameter<bool>("saveMuons");
    saveTaus_         = iConfig.getUntrackedParameter<bool>("saveTaus");
    savePhotons_      = iConfig.getUntrackedParameter<bool>("savePhotons");
    saveMET_          = iConfig.getUntrackedParameter<bool>("saveMET");
    saveGenJets_      = iConfig.getUntrackedParameter<bool>("saveGenJets");
    saveGenParticles_ = iConfig.getUntrackedParameter<bool>("saveGenParticles");

    ecalTPFilterTag_    = iConfig.getUntrackedParameter<edm::InputTag>("ecalTPFilterTag");
    ecalBEFilterTag_    = iConfig.getUntrackedParameter<edm::InputTag>("ecalBEFilterTag");
    hcalHBHEFilterTag_  = iConfig.getUntrackedParameter<edm::InputTag>("hcalHBHEFilterTag");
    hcalLaserFilterTag_ = iConfig.getUntrackedParameter<edm::InputTag>("hcalLaserFilterTag");

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

    int vtxCount, jetCount, jptCount, metCount, muCount, pfMuCount, eleCount, photonCount, pfPhotonCount, tauCount, genCount, genPartCount;
    vtxCount = jetCount = jptCount = metCount = muCount = pfMuCount = eleCount = photonCount = pfPhotonCount = tauCount = genCount = genPartCount = 0;


    /////////////////////////////////////
    // Get PF candidates for later use //
    /////////////////////////////////////


    Handle<PFCandidateCollection> pfCands;
    iEvent.getByLabel(partFlowTag_,pfCands);
    const  PFCandidateCollection thePfColl = *(pfCands.product());


    //////////////////////////
    //Get vertex information//
    //////////////////////////


    Handle<reco::VertexCollection> primaryVtcs;
    iEvent.getByLabel(primaryVtxTag_, primaryVtcs);

    for(VertexCollection::const_iterator iVtx = primaryVtcs->begin(); iVtx!= primaryVtcs->end(); ++iVtx){
        reco::Vertex myVtx = reco::Vertex(*iVtx);
        if(!myVtx.isValid() || myVtx.isFake()) continue;
        TCPrimaryVtx* vtxCon = new ((*primaryVtx)[vtxCount]) TCPrimaryVtx;
        vtxCon->SetXYZ(myVtx.x(), myVtx.y(), myVtx.z());
        vtxCon->SetNDof(myVtx.ndof());
        vtxCon->SetChi2(myVtx.chi2());
        vtxCon->SetNtracks(myVtx.nTracks()); 
        vtxCon->SetSumPt2Trks(sumPtSquared(myVtx));
        vtxCon->SetIsFake(myVtx.isFake());
        ++vtxCount;
    }

    unsigned ivtx = 0;
    VertexRef myVtxRef(primaryVtcs, ivtx);

    ///////////////////////
    //get jet information//
    ///////////////////////

    Handle<double> rhoCorr;
    iEvent.getByLabel(rhoCorrTag_, rhoCorr);
    rhoFactor = (float)(*rhoCorr);

    Handle<double> rho25Corr;
    iEvent.getByLabel(rho25CorrTag_, rho25Corr);
    rho25Factor = (float)(*rho25Corr);

    //cout<<" RHOS. In eta 4.4 = "<<rhoFactor<<"   in eta25 "<<rho25Factor<<endl;

    if(saveJets_){

        Handle<vector<pat::Jet> > jets;
        iEvent.getByLabel(jetTag_, jets);

        for (vector<pat::Jet>::const_iterator iJet = jets->begin(); iJet != jets->end(); ++iJet) {

            if (iJet->pt() < 10.) continue;

            TCJet* jetCon = new ((*recoJets)[jetCount]) TCJet;

            /*
               if (jetCount == 0) {
               cout << "Uncorrected jet pt: " << iJet->correctedJet(0).pt() 
               << ", corrected jet pt: " << iJet->correctedJet(1).pt() 
               << ", corrected jet pt: " << iJet->correctedJet(2).pt() 
               << ", corrected jet pt: " << iJet->correctedJet(3).pt()
               <<endl;
               if(isRealData)	
               cout<< ", corrected jet pt: " << iJet->correctedJet(4).pt() 
               << endl; 
               }
             */
            //
            jetCon->SetPxPyPzE(iJet->px(), iJet->py(), iJet->pz(), iJet->energy());
            jetCon->SetVtx(0., 0., 0.);
            //cout<<"  jetCon object pt = "<<jetCon->Pt()<<endl;
            jetCon->SetChHadFrac(iJet->chargedHadronEnergyFraction());
            jetCon->SetNeuHadFrac(iJet->neutralHadronEnergyFraction());
            jetCon->SetChEmFrac(iJet->chargedEmEnergyFraction());
            jetCon->SetNeuEmFrac(iJet->neutralEmEnergyFraction());
            jetCon->SetNumConstit(iJet->chargedMultiplicity() + iJet->neutralMultiplicity());
            jetCon->SetNumChPart(iJet->chargedMultiplicity());

            jetCon->SetJetFlavor(iJet->partonFlavour());

            jetCon->SetUncertaintyJES(-1);

            jetCon->SetBDiscriminatorMap("TCHE", iJet->bDiscriminator("trackCountingHighEffBJetTags"));
            jetCon->SetBDiscriminatorMap("TCHP", iJet->bDiscriminator("trackCountingHighPurBJetTags"));
            jetCon->SetBDiscriminatorMap("SSVHE", iJet->bDiscriminator("simpleSecondaryVertexHighEffBJetTags"));
            jetCon->SetBDiscriminatorMap("JPB", iJet->bDiscriminator("jetProbabilityBJetTags"));


            /////////////////////////
            // Associate to vertex //
            /////////////////////////

            associateJetToVertex(*iJet, primaryVtcs, jetCon);

            ++jetCount;
        }   

        // For VBF analysis...

        /*
           Handle<reco::JPTJetCollection> jptJets;
           iEvent.getByLabel("ak5JPTJetsL1L2L3", jptJets);

           for (reco::JPTJetCollection::const_iterator iJet = jptJets->begin(); iJet != jptJets->end(); ++iJet) {

        // Perform some pre-cleaning
        edm::RefToBase<reco::Jet> jptjetRef = iJet->getCaloJetRef();
        reco::CaloJet const * rawCaloJet = dynamic_cast<reco::CaloJet const*>(&*jptjetRef);

        if (
        iJet->pt() < 10.
        || (iJet->chargedEmEnergyFraction() + iJet->neutralEmEnergyFraction()) > 0.01
        || rawCaloJet->n90() < 2
        //|| (*jetsID)[(*iJet).getCaloJetRef()].fHPD > 0.98 
        ) continue;

        TCJet* jetCon = new ((*recoJPT)[jptCount]) TCJet;

        //if (jptCount == 0) cout << iJet->pt() << endl;

        jetCon->SetPxPyPzE(iJet->px(), iJet->py(), iJet->pz(), iJet->energy());
        jetCon->SetVtx(0., 0., 0.);

        //jetCon->SetIDMap("Zch", iJet->getSpecific().Zch);
        jetCon->SetChHadFrac(iJet->chargedHadronEnergyFraction());
        jetCon->SetNeuHadFrac(iJet->neutralHadronEnergyFraction());
        jetCon->SetChEmFrac(iJet->chargedEmEnergyFraction());
        jetCon->SetNeuEmFrac(iJet->neutralEmEnergyFraction());
        jetCon->SetNumConstit(iJet->chargedMultiplicity());// + iJet->neutralMultiplicity());
        jetCon->SetNumChPart(iJet->chargedMultiplicity());

        ++jptCount;
        }
         */
    }

    /////////////
    // Get MET //
    /////////////

    if (saveMET_){

        //Handle<vector<pat::MET> > MET;
        //iEvent.getByLabel(metTag_, MET);
        //vector<pat::MET>::const_iterator met = MET->begin();

        Handle<PFMETCollection> MET;
        iEvent.getByLabel(metTag_, MET);
        PFMETCollection::const_iterator met = MET->begin();

        if (MET->begin() != MET->end()) {
            recoMET->SetSumEt(met->sumEt());
            recoMET->SetMagPhi(met->et(), met->phi());

            // PF specififc methods
            recoMET->SetMuonFraction(met->MuonEtFraction());
            recoMET->SetNeutralHadronFraction(met->NeutralHadEtFraction());
            recoMET->SetNeutralEMFraction(met->NeutralEMFraction());
            recoMET->SetChargedHadronFraction(met->ChargedHadEtFraction());
            recoMET->SetChargedEMFraction(met->ChargedEMEtFraction());

        }
    }

    ///////////////
    // Get muons //
    ///////////////

    if (saveMuons_) {

        Handle<vector<reco::Muon> > muons;
        iEvent.getByLabel(muonTag_, muons);

        for (vector<reco::Muon>::const_iterator iMuon = muons->begin(); iMuon != muons->end(); ++iMuon) {
            if (!(iMuon->isGlobalMuon() && iMuon->isTrackerMuon())) continue;

            TCMuon* muCon = new ((*recoMuons)[muCount]) TCMuon;

            muCon->SetPxPyPzE(iMuon->px(), iMuon->py(), iMuon->pz(), iMuon->energy());
            muCon->SetVtx(iMuon->globalTrack()->vx(),iMuon->globalTrack()->vy(),iMuon->globalTrack()->vz());
            muCon->SetPtError(iMuon->globalTrack()->ptError());
            muCon->SetCharge(iMuon->charge());

            // Muon ID variables
            muCon->SetIsPF(iMuon->isPFMuon());
            muCon->SetIsGLB(iMuon->isGlobalMuon());
            muCon->SetIsTRK(iMuon->isTrackerMuon());
            muCon->SetNormalizedChi2(iMuon->globalTrack()->normalizedChi2());
            muCon->SetCaloComp(iMuon->caloCompatibility());
            muCon->SetSegComp(muon::segmentCompatibility(*iMuon));

            muCon->SetNumberOfMatches(iMuon->numberOfMatches());
            muCon->SetNumberOfValidPixelHits(iMuon->globalTrack()->hitPattern().numberOfValidPixelHits());
            muCon->SetNumberOfValidTrackerHits(iMuon->globalTrack()->hitPattern().numberOfValidTrackerHits()); 
            muCon->SetNumberOfValidMuonHits(iMuon->globalTrack()->hitPattern().numberOfValidMuonHits());
            muCon->SetNumberOfLostPixelHits(iMuon->globalTrack()->hitPattern().numberOfLostPixelHits());
            muCon->SetNumberOfLostTrackerHits(iMuon->globalTrack()->hitPattern().numberOfLostTrackerHits());


            // Set isolation map values
            // Detector-based isolation

            muCon->SetIsoMap("NTracks_R03", iMuon->isolationR03().nTracks);
            muCon->SetIsoMap("EmIso_R03", iMuon->isolationR03().emEt);
            muCon->SetIsoMap("HadIso_R03", iMuon->isolationR03().hadEt);
            muCon->SetIsoMap("SumPt_R03", iMuon->isolationR03().sumPt);

            muCon->SetIsoMap("NTracks_R05", iMuon->isolationR05().nTracks);
            muCon->SetIsoMap("EmIso_R05", iMuon->isolationR05().emEt);
            muCon->SetIsoMap("HadIso_R05", iMuon->isolationR05().hadEt);
            muCon->SetIsoMap("SumPt_R05", iMuon->isolationR05().sumPt);

            // PF-based isolation
            muCon->SetIsoMap("pfChargedPt_R03", iMuon->pfIsolationR03().sumChargedParticlePt);
            muCon->SetIsoMap("pfChargedHadronPt_R03", iMuon->pfIsolationR03().sumChargedHadronPt);
            muCon->SetIsoMap("pfPhotonEt_R03", iMuon->pfIsolationR03().sumPhotonEt);
            muCon->SetIsoMap("pfNeutralHadronEt_R03", iMuon->pfIsolationR03().sumNeutralHadronEt);
            muCon->SetIsoMap("pfPUPt_R03", iMuon->pfIsolationR03().sumPUPt);

            muCon->SetIsoMap("pfChargedPt_R04", iMuon->pfIsolationR04().sumChargedParticlePt);
            muCon->SetIsoMap("pfChargedHadronPt_R04", iMuon->pfIsolationR04().sumChargedHadronPt);
            muCon->SetIsoMap("pfPhotonEt_R04", iMuon->pfIsolationR04().sumPhotonEt);
            muCon->SetIsoMap("pfNeutralHadronEt_R04", iMuon->pfIsolationR04().sumNeutralHadronEt);
            muCon->SetIsoMap("pfPUPt_R04", iMuon->pfIsolationR04().sumPUPt);

            muCount++;
        }
    }


    ///////////////////
    // Get electrons //
    ///////////////////


    edm::Handle<reco::ConversionCollection> hConversions;
    iEvent.getByLabel("allConversions", hConversions);

    if (saveElectrons_) {

        Handle<reco::GsfElectronCollection > electrons;
        iEvent.getByLabel(electronTag_, electrons);

        for (vector<reco::GsfElectron>::const_iterator iElectron = electrons->begin(); iElectron != electrons->end(); ++iElectron) {
            if (iElectron->pt() < 10) continue;

            TCElectron* eleCon = new ((*recoElectrons)[eleCount]) TCElectron;

            //cout << fabs(iElectron->eta() - iElectron->superCluster()->eta()) << endl;

            // Basic physics object info
            eleCon->SetPxPyPzE(iElectron->px(), iElectron->py(), iElectron->pz(), iElectron->energy());
            eleCon->SetVtx(iElectron->gsfTrack()->vx(),iElectron->gsfTrack()->vy(),iElectron->gsfTrack()->vz());
            eleCon->SetCharge(iElectron->charge());

            // Fiducial variables
            eleCon->SetIsEB(iElectron->isEB());
            eleCon->SetIsEE(iElectron->isEE());
            eleCon->SetIsInGap(iElectron->isGap());

            // Electron ID variables
            eleCon->SetHadOverEm(iElectron->hadronicOverEm());
            eleCon->SetDphiSuperCluster(iElectron->deltaPhiSuperClusterTrackAtVtx());
            eleCon->SetDetaSuperCluster(iElectron->deltaEtaSuperClusterTrackAtVtx());
            eleCon->SetSigmaIetaIeta(iElectron->sigmaIetaIeta());
            eleCon->SetFBrem(iElectron->fbrem());
            eleCon->SetEOverP(iElectron->eSuperClusterOverP());
            eleCon->SetSCEta(iElectron->superCluster()->eta());

            eleCon->SetPtError(iElectron->gsfTrack()->ptError());
            eleCon->SetNormalizedChi2(iElectron->gsfTrack()->normalizedChi2());

            eleCon->SetNumberOfValidPixelHits(iElectron->gsfTrack()->hitPattern().numberOfValidPixelHits());
            eleCon->SetNumberOfValidTrackerHits(iElectron->gsfTrack()->hitPattern().numberOfValidTrackerHits());
            eleCon->SetNumberOfLostPixelHits(iElectron->gsfTrack()->hitPattern().numberOfLostPixelHits());
            eleCon->SetNumberOfLostTrackerHits(iElectron->gsfTrack()->hitPattern().numberOfLostTrackerHits());


            // Electron Iso variables
            eleCon->SetIsoMap("EmIso_R03", iElectron->dr03EcalRecHitSumEt());
            eleCon->SetIsoMap("HadIso_R03", iElectron->dr03HcalTowerSumEt());
            eleCon->SetIsoMap("SumPt_R03", iElectron->dr03TkSumPt());

            eleCon->SetIsoMap("EmIso_R04", iElectron->dr04EcalRecHitSumEt());
            eleCon->SetIsoMap("HadIso_R04", iElectron->dr04HcalTowerSumEt());
            eleCon->SetIsoMap("SumPt_R04", iElectron->dr04TkSumPt());

            eleCon->SetIsoMap("pfPhotonEt_R03", iElectron->pfIsolationVariables().photonIso);
            eleCon->SetIsoMap("pfChargedHadron_R03", iElectron->pfIsolationVariables().chargedHadronIso);
            eleCon->SetIsoMap("pfNeutralHadron_R03", iElectron->pfIsolationVariables().neutralHadronIso);

            eleIsolator.fGetIsolation(&(*iElectron), &thePfColl, myVtxRef, primaryVtcs);
            eleCon->SetIsoMap("pfChIso_R04",eleIsolator.getIsolationCharged());
            eleCon->SetIsoMap("pfNeuIso_R04",eleIsolator.getIsolationNeutral());
            eleCon->SetIsoMap("pfPhoIso_R04",eleIsolator.getIsolationPhoton());

            // Effective area for rho PU corrections (not sure if needed)
            float AEff03 = ElectronEffectiveArea::GetElectronEffectiveArea(ElectronEffectiveArea::kEleGammaAndNeutralHadronIso03, iElectron->eta(), ElectronEffectiveArea::kEleEAData2012);
            float AEff04 = ElectronEffectiveArea::GetElectronEffectiveArea(ElectronEffectiveArea::kEleGammaAndNeutralHadronIso04, iElectron->eta(), ElectronEffectiveArea::kEleEAData2012);
            eleCon->SetIdMap("EffArea_R03", AEff03);
            eleCon->SetIdMap("EffArea_R04", AEff04);

            // Conversion information
            bool convVeto = !(ConversionTools::hasMatchedConversion(*iElectron,hConversions,vertexBeamSpot.position()));
            eleCon->SetConversionVeto(convVeto);
            eleCon->SetConversionMissHits(iElectron->gsfTrack()->trackerExpectedHitsInner().numberOfHits());

            // EID maps for VBTF working points -- probably not needed anymore
            //eleCon->SetCutLevel(iElectron->electronID("eidVBTF95"), 95);
            //eleCon->SetCutLevel(iElectron->electronID("eidVBTF90"), 90);
            //eleCon->SetCutLevel(iElectron->electronID("eidVBTF85"), 85);
            //eleCon->SetCutLevel(iElectron->electronID("eidVBTF80"), 80);
            //eleCon->SetCutLevel(iElectron->electronID("eidVBTF70"), 70);
            //eleCon->SetCutLevel(iElectron->electronID("eidVBTF60"), 60);

            // Add electron MVA ID and ISO when done -- needs work
            //electronMVA(vtxCollection, iElectron);

            eleCount++;
        }
    }


    /////////////////
    // Get photons //
    /////////////////


    if (savePhotons_) {

        Handle<vector<reco::Photon> > photons;
        iEvent.getByLabel(photonTag_, photons);

        edm::Handle<reco::GsfElectronCollection> hElectrons;
        iEvent.getByLabel("gsfElectrons", hElectrons);

        for (vector<reco::Photon>::const_iterator iPhoton = photons->begin(); iPhoton != photons->end() ; ++iPhoton) {

            TCPhoton* myPhoton = new ((*recoPhotons)[photonCount]) TCPhoton();
            myPhoton->SetPxPyPzE(iPhoton->px(), iPhoton->py(), iPhoton->pz(), iPhoton->p());
            myPhoton->SetVtx(iPhoton->vx(), iPhoton->vy(), iPhoton->vz());

            // ID variables
            myPhoton->SetHadOverEm(iPhoton->hadronicOverEm());
            myPhoton->SetSigmaIEtaIEta(iPhoton->sigmaIetaIeta());
            myPhoton->SetR9(iPhoton->r9());
            myPhoton->SetTrackVeto(iPhoton->hasPixelSeed());

            myPhoton->SetEtaSC(iPhoton->superCluster()->eta());
            myPhoton->SetEnergySC(iPhoton->superCluster()->energy());

            // detector-based isolation
            myPhoton->SetIsoMap("EmIso_R03", (iPhoton->ecalRecHitSumEtConeDR03()));
            myPhoton->SetIsoMap("HadIso_R03", (iPhoton->hcalTowerSumEtConeDR03()));
            myPhoton->SetIsoMap("TrkIso_R03", (iPhoton->trkSumPtHollowConeDR03()));

            myPhoton->SetIsoMap("EmIso_R04", (iPhoton->ecalRecHitSumEtConeDR04()));
            myPhoton->SetIsoMap("HadIso_R04", (iPhoton->hcalTowerSumEtConeDR04()));
            myPhoton->SetIsoMap("TrkIso_R04", (iPhoton->trkSumPtHollowConeDR04()));

            // PF Iso for photons
            phoIsolator.fGetIsolation(&(*iPhoton),&thePfColl, myVtxRef, primaryVtcs);
            myPhoton->SetIsoMap("chIso03",phoIsolator.getIsolationCharged());
            myPhoton->SetIsoMap("nhIso03",phoIsolator.getIsolationNeutral());
            myPhoton->SetIsoMap("phIso03",phoIsolator.getIsolationPhoton());

            // Hcal isolation for 2012
            //myPhoton->SetIsoMap("HadIso_R03",iPhoton->hcalTowerSumEtConeDR03() + 
            //        (iPhoton->hadronicOverEm() - iPhoton->hadTowOverEm())*iPhoton->superCluster()->energy()/cosh(iPhoton->superCluster()->eta()));
            //myPhoton->SetIsoMap("HadIso_R04",iPhoton->hcalTowerSumEtConeDR04() + 
            //        (iPhoton->hadronicOverEm() - iPhoton->hadTowOverEm())*iPhoton->superCluster()->energy()/cosh(iPhoton->superCluster()->eta()));


            //Conversion info
            bool passElectronVeto = !(ConversionTools::hasMatchedPromptElectron(iPhoton->superCluster(), hElectrons, hConversions, vertexBeamSpot.position()));
            myPhoton->SetConversionVeto(passElectronVeto);

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

            if (!iTau->isPFTau() 
                    or iTau->signalPFChargedHadrCands().size() < 1 
                    or iTau->pt() < 10
               ) continue;

            TCTau* tauCon = new ((*recoTaus)[tauCount]) TCTau;

            tauCon->SetNChHad(iTau->signalPFChargedHadrCands().size());
            tauCon->SetNGamma (iTau->signalPFGammaCands().size());
            tauCon->SetNNeutrHad (iTau->signalPFNeutrHadrCands().size());
            tauCon->SetCharge(iTau->charge());
            tauCon->SetDecayMode(iTau->decayMode());

            tauCon->SetPxPyPzE(iTau->px(),iTau->py(),iTau->pz(),iTau->energy());
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
        } 
    } 


    ////////////////////////
    // Get gen-level info //
    ////////////////////////


    if (!isRealData) {

        Handle<GenEventInfoProduct> GenEventInfoHandle;
        iEvent.getByLabel("generator", GenEventInfoHandle);

        evtWeight = ptHat = qScale = -1;

        if (GenEventInfoHandle.isValid()) {
            //qScale       = GenEventInfoHandle->qScale();
            //evtWeight    = GenEventInfoHandle->weight();
            ptHat        = (GenEventInfoHandle->hasBinningValues() ? GenEventInfoHandle->binningValues()[0] : 0.0);
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

                //cout<<iPV->getTrueNumInteractions()<<endl;
            }
        }

        //////////////////////
        // Get genParticles //
        //////////////////////

        if (saveGenParticles_) {
            Handle<GenParticleCollection> genParticleColl;
            iEvent.getByLabel("genParticles", genParticleColl);

            for (GenParticleCollection::const_iterator iGenPart = genParticleColl->begin(); iGenPart != genParticleColl->end(); ++iGenPart) {
                const reco::GenParticle myParticle = reco::GenParticle(*iGenPart);

                ////  Leptons and photons and b's, (oh my)
                if (
                        myParticle.pt() > 8 
                        && (
                            (abs(myParticle.pdgId()) >= 11 && abs(myParticle.pdgId()) <= 16) 
                            || myParticle.pdgId() == 22 
                            || abs(myParticle.pdgId()) == 5 
                           )
                   ) {

                    TCGenParticle* genCon = new ((*genParticles)[genPartCount]) TCGenParticle;
                    genCon->SetPxPyPzE(myParticle.px(), myParticle.py(), myParticle.pz(), myParticle.energy() );
                    genCon->SetVtx(myParticle.vx(), myParticle.vy(), myParticle.vz());
                    genCon->SetCharge(myParticle.charge());
                    genCon->SetPDGId(myParticle.pdgId());
                    genCon->SetMother(myParticle.mother()->pdgId());
                    genCon->SetStatus(myParticle.status());
                    if (myParticle.mother()->numberOfMothers() != 0) genCon->SetGrandmother(myParticle.mother()->mother()->pdgId());
                    ++genPartCount;
                }

                //// Z's, W's, H's, and now big juicy Gravitons
                if (
                        abs(myParticle.pdgId()) == 23 
                        || abs(myParticle.pdgId()) == 24 
                        || abs(myParticle.pdgId()) == 25 
                        || abs(myParticle.pdgId()) == 39
                   ){


                    TCGenParticle* genCon = new ((*genParticles)[genPartCount]) TCGenParticle;
                    genCon->SetPxPyPzE(myParticle.px(), myParticle.py(), myParticle.pz(), myParticle.energy() );
                    genCon->SetVtx(myParticle.vx(), myParticle.vy(), myParticle.vz() );
                    genCon->SetCharge(myParticle.charge());
                    genCon->SetPDGId(myParticle.pdgId());
                    genCon->SetMother(myParticle.mother()->pdgId());
                    genCon->SetStatus(myParticle.status());
                    ++genPartCount;
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
                    jetCon->SetPxPyPzE(myJet.px(), myJet.py(), myJet.pz(), myJet.energy());
                    jetCon->SetHadEnergy(myJet.hadEnergy());
                    jetCon->SetEmEnergy(myJet.emEnergy());
                    jetCon->SetInvEnergy(myJet.invisibleEnergy());
                    jetCon->SetAuxEnergy(myJet.auxiliaryEnergy());
                    jetCon->SetNumConstit(myJet.getGenConstituents().size());
                    jetCon->SetJetFlavor(0);
                }
                ++genCount;	
            }
        }
    }


    ///////////////////
    // Noise filters //
    ///////////////////

    //if (isRealData) {

    myNoiseFilters.isScraping = isFilteredOutScraping(iEvent, iSetup, 10, 0.25);

    Handle<bool> hcalNoiseFilterHandle;
    iEvent.getByLabel(hcalHBHEFilterTag_, hcalNoiseFilterHandle);
    if (hcalNoiseFilterHandle.isValid())  myNoiseFilters.isNoiseHcalHBHE = !(Bool_t)(*hcalNoiseFilterHandle);
    else LogWarning("Filters")<<"hcal noise NOT valid  ";

    Handle<bool> hcalLaserFilterHandle;
    iEvent.getByLabel(hcalLaserFilterTag_, hcalLaserFilterHandle);
    if (hcalLaserFilterHandle.isValid())  myNoiseFilters.isNoiseHcalLaser = !(Bool_t)(*hcalLaserFilterHandle);
    else LogWarning("Filters")<<"hcal Laser NOT valid  ";

    Handle<bool> ecalTPFilterHandle;
    iEvent.getByLabel(ecalTPFilterTag_, ecalTPFilterHandle);
    if (ecalTPFilterHandle.isValid())  myNoiseFilters.isNoiseEcalTP = !(Bool_t)(*ecalTPFilterHandle);
    else LogWarning("Filters")<<"Ecal TP NOT valid  ";

    Handle<bool> ecalBEFilterHandle;
    iEvent.getByLabel(ecalBEFilterTag_, ecalBEFilterHandle);
    if (ecalBEFilterHandle.isValid())  myNoiseFilters.isNoiseEcalBE = !(Bool_t)(*ecalBEFilterHandle);
    else LogWarning("Filters")<<"Ecal BE NOT valid  ";

    edm::Handle<BeamHaloSummary> TheBeamHaloSummary;
    iEvent.getByLabel("BeamHaloSummary",TheBeamHaloSummary);
    const BeamHaloSummary TheSummary = (*TheBeamHaloSummary.product() );
    if (!TheBeamHaloSummary.isValid()) LogWarning("Filters")<<"The Summary (for CSC halo) NOT valid  ";

    myNoiseFilters.isCSCTightHalo = TheSummary.CSCTightHaloId();
    myNoiseFilters.isCSCLooseHalo = TheSummary.CSCLooseHaloId();


    //LogWarning("Filters")<<"\n csc1  "<< myNoiseFilters.isCSCTightHalo<<"  csc2  "<<myNoiseFilters.isCSCLooseHalo
    //  <<" isNoiseHcal HBHE "<<myNoiseFilters.isNoiseHcalHBHE<<"  laser "<<myNoiseFilters.isNoiseHcalLaser<<"\n"
    //  <<" ecal TP  "<<myNoiseFilters.isNoiseEcalTP<<"   ecal BE  "<<myNoiseFilters.isNoiseEcalBE;

    //}

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

    edm::Handle<trigger::TriggerEvent> triggerEvents;
    iEvent.getByLabel("hltTriggerSummaryAOD",triggerEvents);
    trigger::TriggerObjectCollection triggerObjCol = triggerEvents->getObjects();
    int triggerCount = 0;

    for(trigger::TriggerObjectCollection::const_iterator iTrigObj = triggerObjCol.begin(); iTrigObj != triggerObjCol.end(); ++iTrigObj) { 
        TCTriggerObject * thisTrig = new ((*triggerObjects)[triggerCount]) TCTriggerObject;
        thisTrig->SetPxPyPzE(iTrigObj->px(), iTrigObj->py(), iTrigObj->pz(), iTrigObj->energy());
        thisTrig->SetId(iTrigObj->id());
        ++triggerCount;
    }

    ++nEvents;

    if (eleCount > 0 || muCount > 0) eventTree -> Fill(); // possibly specify a cut in configuration

    primaryVtx    -> Clear("C");
    recoJets      -> Clear("C");
    recoJPT       -> Clear("C");
    recoMuons     -> Clear("C");
    recoElectrons -> Clear("C");
    recoTaus      -> Clear("C");
    recoPhotons   -> Clear("C");
    //pfPhotons   -> Clear("C");
    triggerObjects-> Clear("C");
    genJets       -> Clear("C");
    genParticles  -> Clear("C");
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
    //pfPhotons      = new TClonesArray("TCPhoton");
    triggerObjects = new TClonesArray("TCTriggerObject");
    genJets        = new TClonesArray("TCGenJet");
    genParticles   = new TClonesArray("TCGenParticle");
    beamSpot       = new TVector3();
    recoMET        = 0;

    eventTree->Branch("recoJets",&recoJets, 6400, 0);
    eventTree->Branch("recoJPT",&recoJPT, 6400, 0);
    eventTree->Branch("recoElectrons",&recoElectrons, 6400, 0);
    eventTree->Branch("recoMuons",&recoMuons, 6400, 0);
    eventTree->Branch("recoTaus",&recoTaus, 6400, 0);
    eventTree->Branch("recoPhotons",&recoPhotons, 6400, 0);
    //eventTree->Branch("pfPhotons",&pfPhotons, 6400, 0);
    eventTree->Branch("recoMET", &recoMET, 6400, 0);
    eventTree->Branch("triggerObjects", &triggerObjects, 6400, 0);
    eventTree->Branch("genJets",&genJets, 6400, 0);
    eventTree->Branch("genParticles",&genParticles, 6400, 0);

    eventTree->Branch("primaryVtx",&primaryVtx, 6400, 0);
    eventTree->Branch("beamSpot", &beamSpot, 6400, 0);
    eventTree->Branch("nPUVertices", &nPUVertices, "nPUVertices/I");
    eventTree->Branch("nPUVerticesTrue", &nPUVerticesTrue, "nPUVerticesTrue/F");

    eventTree->Branch("isRealData",&isRealData, "isRealData/O");
    eventTree->Branch("runNumber",&runNumber, "runNumber/i");
    eventTree->Branch("eventNumber",&eventNumber, "eventNumber/l");
    eventTree->Branch("lumiSection",&lumiSection, "lumiSection/i");
    eventTree->Branch("bunchCross",&bunchCross, "bunchCross/i");

    eventTree->Branch("ptHat",&ptHat, "ptHat/F");
    eventTree->Branch("qScale", &qScale, "qScale/F");
    eventTree->Branch("evtWeight", &evtWeight, "evtWeight/F");
    eventTree->Branch("rhoFactor",&rhoFactor, "rhoFactor/F");
    eventTree->Branch("rho25Factor",&rho25Factor, "rho25Factor/F");
    eventTree->Branch("triggerStatus",&triggerStatus, "triggerStatus/l");
    eventTree->Branch("hltPrescale",hltPrescale, "hltPrescale[64]/i");

    eventTree->Branch("NoiseFilters", &myNoiseFilters.isScraping, "isScraping/O:isNoiseHcalHBHE:isNoiseHcalLaser:isNoiseEcalTP:isNoiseEcalBE:isCSCTightHalo:isCSCLooseHalo");

    runTree->Branch("deliveredLumi",&deliveredLumi, "deliveredLumi/F");
    runTree->Branch("recordedLumi",&recordedLumi, "recordedLumi/F");
    runTree->Branch("runNumber",&runNumber, "runNumber/i");

    jobTree->Branch("savedTriggerNames",savedTriggerNames, "savedTriggerNames[64]/C");
    jobTree->Branch("nEvents",&nEvents, "nEvents/i");

    // Initialize HLT prescales //

    for (int i = 0; i < (int)(sizeof(hltPrescale)/sizeof(int)); ++i) hltPrescale[i] = 1;

    // Start counting number of events per job //
    nEvents = 0;

    // Photon Iso maker init
    phoIsolator.initializePhotonIsolation(kTRUE);
    // Below should be default from initializePhotonIsolation, but we can set them to be sure
    phoIsolator.setConeSize(0.3);
    phoIsolator.setDeltaRVetoBarrelCharged(0.02);
    phoIsolator.setDeltaRVetoEndcapCharged(0.02);
    phoIsolator.setRectangleDeltaEtaVetoBarrelPhotons(0.015);
    phoIsolator.setDeltaRVetoEndcapPhotons(0.07);

    // Initialize Electron MVA nonsense
    eleIsolator.initializeElectronIsolation(kTRUE);
    eleIsolator.setConeSize(0.4);
                                    
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


bool ntupleProducer::associateJetToVertex(pat::Jet inJet, Handle<reco::VertexCollection> vtxCollection, TCJet *outJet)
{
    if(fabs(inJet.eta()) > 2.5){
        outJet->SetVtxSumPtFrac(-1);
        outJet->SetVtxSumPt(-1);
        outJet->SetVtxTrackFrac(-1);
        outJet->SetVtxNTracks(-1);
        outJet->SetVtxSumPtIndex(0);
        outJet->SetVtxCountIndex(0);

        return false;
    }

    vector<float>  associatedTrackSumPt;
    vector<float>  associatedTrackCount;
    vector<const reco::Track*> jetTracks;
    float sumTrackX, sumTrackY, sumTrackZ, sumTrackPt;
    int   nJetTracks = 0;
    int   vCount = 0;

    sumTrackX = sumTrackY = sumTrackZ  = sumTrackPt = 0;

    const reco::TrackRefVector &tracks = inJet.associatedTracks(); 
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

    return true;
}

bool ntupleProducer::electronMVA(Handle<reco::VertexCollection> vtxCollection, vector<pat::Electron>::const_iterator iElectron)
{
    if (vtxCollection->size() != 0) return false;

    /*
       edm::ESHandle<TransientTrackBuilder> builder;
       iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", builder);
       TransientTrackBuilder thebuilder = *(builder.product());

       InputTag  reducedEBRecHitCollection(string("reducedEcalRecHitsEB"));
       InputTag  reducedEERecHitCollection(string("reducedEcalRecHitsEE"));

       EcalClusterLazyTools lazyTools(iEvent, iSetup, reducedEBRecHitCollection, reducedEERecHitCollection);

       pv = &*vtxCollection->begin(); 
       double myMVANonTrigMethod1 = myMVANonTrig->mvaValue(*iElectron,*pv,thebuilder,lazyTools,false);
       double myMVATrigMethod1 = myMVATrig->mvaValue(*iElectron,*pv,thebuilder,lazyTools,false);

    //eleCon->SetIdMap("MVATrigMethod1",myMVATrigMethod1);
    //eleCon->SetIdMap("MVANonTrigMethod1",myMVANonTrigMethod1);

    //ID'd electrons to feed into MVAIso

    InputTag gsfEleLabel(string("gsfElectrons"));
    Handle<reco::GsfElectronCollection> theEGammaCollection;
    iEvent.getByLabel(gsfEleLabel,theEGammaCollection);

    reco::GsfElectronCollection identifiedElectrons;
    for (reco::GsfElectronCollection::const_iterator iE = theEGammaCollection->begin(); iE != theEGammaCollection->end(); ++iE) {

    double electronTrackZ = 0;
    if (iE->gsfTrack().isNonnull()) {
    electronTrackZ = iE->gsfTrack()->dz(pv->position());
    } else if (iE->closestCtfTrackRef().isNonnull()) {
    electronTrackZ = iE->closestCtfTrackRef()->dz(pv->position());
    }    
    if(fabs(electronTrackZ) > 0.2)  continue;


    if(fabs(iE->superCluster()->eta())<1.479) {     
    if(iE->pt() > 20) {
    if(
    iE->sigmaIetaIeta()       > 0.01
    || fabs(iE->deltaEtaSuperClusterTrackAtVtx()) > 0.007
    || fabs(iE->deltaPhiSuperClusterTrackAtVtx()) > 0.8
    || iE->hadronicOverEm()       > 0.15
    )  continue;    
    } else {
    if(iE->sigmaIetaIeta()       > 0.012)  continue;
    if(fabs(iE->deltaEtaSuperClusterTrackAtVtx()) > 0.007) continue;
    if(fabs(iE->deltaPhiSuperClusterTrackAtVtx()) > 0.8)  continue;
    if(iE->hadronicOverEm()       > 0.15) continue;    
    } 
    } else {     
    if(iE->pt() > 20) {
    if(iE->sigmaIetaIeta()       > 0.03)  continue;
    if(fabs(iE->deltaEtaSuperClusterTrackAtVtx()) > 0.010) continue;
    if(fabs(iE->deltaPhiSuperClusterTrackAtVtx()) > 0.8)  continue;
    } else {
    if(iE->sigmaIetaIeta()       > 0.032)  continue;
    if(fabs(iE->deltaEtaSuperClusterTrackAtVtx()) > 0.010) continue;
    if(fabs(iE->deltaPhiSuperClusterTrackAtVtx()) > 0.8)  continue;
    }
    }
    identifiedElectrons.push_back(*iE);
    }

    Handle<reco::MuonCollection> hMuonProduct;
    iEvent.getByLabel("muons", hMuonProduct);  
    const reco::MuonCollection inMuons = *(hMuonProduct.product()); 

    reco::MuonCollection identifiedMuons;
    for (reco::MuonCollection::const_iterator iMuon = inMuons.begin(); iMuon != inMuons.end(); ++iMuon) {
    if (
    iMuon->innerTrack().isNonnull()
    && iMuon->isGlobalMuon() 
        && iMuon->isTrackerMuon();
    && iMuon->innerTrack()->numberOfValidHits() > 11 
        ) identifiedMuons.push_back(*iMuon);
}

Handle<PFCandidateCollection> pfCands;
iEvent.getByLabel("particleFlow", pfCands);
const  PFCandidateCollection pfCanIso = *(pfCands.product());
double isomva = fElectronIsoMVA->mvaValue( *iElectron, *pv, pfCanIso, rhoFactor, ElectronEffectiveArea::kEleEAData2011, identifiedElectrons, identifiedMuons);
*/

return true;
}

//define this as a plug-in
DEFINE_FWK_MODULE(ntupleProducer);

#include "../interface/ntupleProducer.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"

ntupleProducer::ntupleProducer(const ParameterSet& iConfig)
{
    jetTag_           = iConfig.getUntrackedParameter<InputTag>("JetTag");
    metTag_           = iConfig.getUntrackedParameter<InputTag>("METTag");
    muonTag_          = iConfig.getUntrackedParameter<InputTag>("MuonTag");
    electronTag_      = iConfig.getUntrackedParameter<InputTag>("ElectronTag");
    photonTag_        = iConfig.getUntrackedParameter<InputTag>("PhotonTag");
    genJetTag_        = iConfig.getUntrackedParameter<InputTag>("GenJetTag");
    primaryVtxTag_    = iConfig.getUntrackedParameter<InputTag>("PrimaryVtxTag");

    rhoCorrTag_       = iConfig.getUntrackedParameter<InputTag>("rhoCorrTag");
    rho25CorrTag_     = iConfig.getUntrackedParameter<InputTag>("rho25CorrTag");
    rhoMuCorrTag_     = iConfig.getUntrackedParameter<InputTag>("rhoMuCorrTag");

    hlTriggerResults_ = iConfig.getUntrackedParameter<string>("HLTriggerResults","TriggerResults");
    hltProcess_       = iConfig.getUntrackedParameter<string>("hltName");
    triggerPaths_     = iConfig.getUntrackedParameter<vector<string> >("triggers");

    saveJets_         = iConfig.getUntrackedParameter<bool>("saveJets");
    saveElectrons_    = iConfig.getUntrackedParameter<bool>("saveElectrons");
    saveMuons_        = iConfig.getUntrackedParameter<bool>("saveMuons");
    savePhotons_      = iConfig.getUntrackedParameter<bool>("savePhotons");
    saveMET_          = iConfig.getUntrackedParameter<bool>("saveMET");

    saveGenJets_      = iConfig.getUntrackedParameter<bool>("saveGenJets");
    saveGenParticles_ = iConfig.getUntrackedParameter<bool>("saveGenParticles");

    verboseTrigs       = iConfig.getUntrackedParameter<bool>("verboseTrigs");
    verboseMVAs        = iConfig.getUntrackedParameter<bool>("verboseMVAs");

    hcalHBHEFilterTag_  = iConfig.getUntrackedParameter<InputTag>("hcalHBHEFilterTag");
    photonIsoCalcTag_   = iConfig.getParameter<ParameterSet>("photonIsoCalcTag");
}

ntupleProducer::~ntupleProducer()
{

}

//
// member functions
//

// ------------ method called to for each event  ------------
void ntupleProducer::analyze(const Event& iEvent, const EventSetup& iSetup)
{
    // this is for CVS check test
    eventNumber  = iEvent.id().event(); 
    runNumber    = iEvent.id().run();
    lumiSection  = (unsigned int)iEvent.getLuminosityBlock().luminosityBlock();
    bunchCross   = (unsigned int)iEvent.bunchCrossing();
    isRealData   = iEvent.isRealData();

    Handle<BeamSpot> beamSpotHandle;
    iEvent.getByLabel("offlineBeamSpot", beamSpotHandle);
    BeamSpot vertexBeamSpot = *beamSpotHandle;

    beamSpot->SetXYZ(vertexBeamSpot.x0(), vertexBeamSpot.y0(), vertexBeamSpot.z0());

    int vtxCount, jetCount, jptCount, metCount, muCount, pfMuCount, eleCount, photonCount, pfPhotonCount, genCount, genPartCount, trigCount;
    vtxCount = jetCount = jptCount = metCount = muCount = pfMuCount = eleCount = photonCount = pfPhotonCount = genCount = genPartCount = trigCount = 0;


    /////////////////////////////////////
    // Get PF candidates for later use //
    /////////////////////////////////////


    Handle<PFCandidateCollection> pfCands;
    iEvent.getByLabel("particleFlow",pfCands);
    const  PFCandidateCollection thePfColl = *(pfCands.product());

    Handle<PFCandidateCollection> pfCandsEleIso;
    iEvent.getByLabel("pfNoPileUp",pfCandsEleIso);
    const  PFCandidateCollection thePfCollEleIso = *(pfCandsEleIso.product());


    //////////////////////////
    //Get vertex information//
    //////////////////////////


    Handle<VertexCollection> primaryVtcs;
    iEvent.getByLabel(primaryVtxTag_, primaryVtcs);

    for(VertexCollection::const_iterator iVtx = primaryVtcs->begin(); iVtx!= primaryVtcs->end(); ++iVtx){
        Vertex myVtx = Vertex(*iVtx);
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

    Handle<double> rhoMuCorr;
    iEvent.getByLabel(rhoMuCorrTag_, rhoMuCorr);
    rhoMuFactor = (float)(*rhoMuCorr);

    //cout<<" RHOS. In eta 4.4 = "<<rhoFactor<<"   in eta25 "<<rho25Factor<<"  MUs: "<<rhoMuFactor<<endl;

    if(saveJets_){

        // b-taggers //
        Handle<JetTagCollection> bTagCollectionTCHE;
        iEvent.getByLabel("trackCountingHighEffBJetTags", bTagCollectionTCHE);
        const JetTagCollection & bTagsTCHE = *(bTagCollectionTCHE.product());

        Handle<JetTagCollection> bTagCollectionTCHP;
        iEvent.getByLabel("trackCountingHighPurBJetTags", bTagCollectionTCHP);
        const JetTagCollection & bTagsTCHP = *(bTagCollectionTCHP.product()); 

        Handle<JetTagCollection> bTagCollectionSSVHE;
        iEvent.getByLabel("simpleSecondaryVertexHighEffBJetTags", bTagCollectionSSVHE);
        const JetTagCollection & bTagsSSVHE = *(bTagCollectionSSVHE.product());

        Handle<JetTagCollection> bTagCollectionJBP;
        iEvent.getByLabel("jetBProbabilityBJetTags", bTagCollectionJBP);
        const JetTagCollection & bTagsJBP = *(bTagCollectionJBP.product());

        Handle<JetTagCollection> bTagCollectionCSV;
        iEvent.getByLabel("combinedSecondaryVertexBJetTags", bTagCollectionCSV);
        const JetTagCollection & bTagsCSV = *(bTagCollectionCSV.product());

        // Jet flavor matching //
        Handle<JetFlavourMatchingCollection> bJetFlavourMC;
        if (!isRealData) {
            iEvent.getByLabel("JetbyValAlgo", bJetFlavourMC);
        }

        Handle<vector<PFJet> > jets;
        iEvent.getByLabel(jetTag_, jets);

        for (vector<PFJet>::const_iterator iJet = jets->begin(); iJet != jets->end(); ++iJet) {

            if (iJet->pt() < 15.) continue;

            TCJet* jetCon = new ((*recoJets)[jetCount]) TCJet;

            jetCon->SetPxPyPzE(iJet->px(), iJet->py(), iJet->pz(), iJet->energy());
            jetCon->SetVtx(0., 0., 0.);
            jetCon->SetChHadFrac(iJet->chargedHadronEnergyFraction());
            jetCon->SetNeuHadFrac(iJet->neutralHadronEnergyFraction());
            jetCon->SetChEmFrac(iJet->chargedEmEnergyFraction());
            jetCon->SetNeuEmFrac(iJet->neutralEmEnergyFraction());
            jetCon->SetNumConstit(iJet->chargedMultiplicity() + iJet->neutralMultiplicity());
            jetCon->SetNumChPart(iJet->chargedMultiplicity());

            if (!isRealData) {
                int jetFlavor = 0;
                for (JetFlavourMatchingCollection::const_iterator iFlavor = bJetFlavourMC->begin(); iFlavor != bJetFlavourMC->end(); iFlavor++) {
                    if (iFlavor->first.get()->pt() > 10 && iFlavor->first.get()->pt() == iJet->pt()) {
                        jetFlavor = iFlavor->second.getFlavour();
                    }
                }
                //cout << iJet->pt() << "\t" << jetFlavor << endl;

                jetCon->SetJetFlavor(jetFlavor);
            } else {
                jetCon->SetJetFlavor(0);
            }


            jetCon->SetUncertaintyJES(-1);

            jetCon->SetBDiscriminatorMap("TCHE", MatchBTagsToJets(bTagsTCHE, *iJet));
            jetCon->SetBDiscriminatorMap("TCHP", MatchBTagsToJets(bTagsTCHP, *iJet));
            jetCon->SetBDiscriminatorMap("SSVHE", MatchBTagsToJets(bTagsSSVHE, *iJet));
            jetCon->SetBDiscriminatorMap("JBP", MatchBTagsToJets(bTagsJBP, *iJet));
            jetCon->SetBDiscriminatorMap("CSV", MatchBTagsToJets(bTagsCSV, *iJet));


            /////////////////////////
            // Associate to vertex //
            /////////////////////////

            associateJetToVertex(*iJet, primaryVtcs, jetCon);

            ++jetCount;
        }   
    }


    /////////////
    // Get MET //
    /////////////

    if (saveMET_) {

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

        recoMET->SetSignificance(MET->front().significance());
        recoMET->SetSigmaX2(MET->front().getSignificanceMatrix()(0,0));

        //Handle<View<PFMET>> pfMEThandle;
        //iEvent.getByLabel("pfMet", pfMEThandle);
        //recoMET->SetSignificance(pfMEThandle->front().significance());
        //recoMET->SetSigmaX2(pfMEThandle->front().getSignificanceMatrix()(0,0));


        ///////////////
        // Get T0MET //
        ///////////////


        Handle<PFMETCollection> t0MET;
        iEvent.getByLabel("pfType1CorrectedMetType0", t0MET);
        PFMETCollection::const_iterator t0met = t0MET->begin();

        if (t0MET->begin() != t0MET->end()) {
            T0MET->SetSumEt(t0met->sumEt());
            T0MET->SetMagPhi(t0met->et(), t0met->phi());

            // PF specififc methods
            T0MET->SetMuonFraction(t0met->MuonEtFraction());
            T0MET->SetNeutralHadronFraction(t0met->NeutralHadEtFraction());
            T0MET->SetNeutralEMFraction(t0met->NeutralEMFraction());
            T0MET->SetChargedHadronFraction(t0met->ChargedHadEtFraction());
            T0MET->SetChargedEMFraction(t0met->ChargedEMEtFraction());

            //Significance
            float significance = (t0MET->front()).significance();
            float sigmaX2 = (t0MET->front()).getSignificanceMatrix()(0,0);
            T0MET->SetSignificance( significance );
            T0MET->SetSigmaX2( sigmaX2 );

        }

        ///////////////
        // Get T2MET //
        ///////////////


        Handle<PFMETCollection> t2MET;
        iEvent.getByLabel("pfType1p2CorrectedMet", t2MET);
        PFMETCollection::const_iterator t2met = t2MET->begin();

        if (t2MET->begin() != t2MET->end()) {
            T2MET->SetSumEt(t2met->sumEt());
            T2MET->SetMagPhi(t2met->et(), t2met->phi());

            // PF specififc methods
            T2MET->SetMuonFraction(t2met->MuonEtFraction());
            T2MET->SetNeutralHadronFraction(t2met->NeutralHadEtFraction());
            T2MET->SetNeutralEMFraction(t2met->NeutralEMFraction());
            T2MET->SetChargedHadronFraction(t2met->ChargedHadEtFraction());
            T2MET->SetChargedEMFraction(t2met->ChargedEMEtFraction());

        }

        //////////////////
        // Get TrackMET //  
        //////////////////


        Handle<METCollection> trkMET;
        iEvent.getByLabel("trackMet", trkMET);
        METCollection::const_iterator trkmet = trkMET->begin();

        if (trkMET->begin() != trkMET->end()) {
            track_MET->SetSumEt(trkmet->sumEt());
            track_MET->SetMagPhi(trkmet->et(), trkmet->phi());
        }

        ////////////////
        // Get MVAMET // 
        //////////////// 

        Handle<vector<PFMET> > mvaMET;
        iEvent.getByLabel("pfMEtMVA", mvaMET);
        vector<PFMET>::const_iterator mvamet = mvaMET->begin();

        if (mvaMET->begin() != mvaMET->end()) {
            //std::cout << "sumET: " << mvamet->sumEt()<<std::endl;
            //std::cout << "et: " << mvamet->et() << std::endl;
            //std::cout << "phi: " << mvamet->phi() << std::endl;
            mva_MET->SetSumEt(mvamet->sumEt());
            mva_MET->SetMagPhi(mvamet->et(), mvamet->phi());
        }
    }



    ///////////////
    // Get muons //
    ///////////////


    if (saveMuons_) {

        Handle<vector<Muon> > muons;
        iEvent.getByLabel(muonTag_, muons);

        for (vector<Muon>::const_iterator iMuon = muons->begin(); iMuon != muons->end(); ++iMuon) {
            if (!iMuon->isGlobalMuon() || iMuon->pt() < 3.) continue;
            //if (!iMuon->isGlobalMuon()) continue;

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

            muCon->SetTrackLayersWithMeasurement(iMuon->track()->hitPattern().trackerLayersWithMeasurement());
            muCon->SetNumberOfMatchedStations(iMuon->numberOfMatchedStations());
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


    Handle<ConversionCollection> hConversions;
    iEvent.getByLabel("allConversions", hConversions);

    if (saveElectrons_) {

        Handle<GsfElectronCollection > electrons;
        iEvent.getByLabel(electronTag_, electrons);

        for (vector<GsfElectron>::const_iterator iElectron = electrons->begin(); iElectron != electrons->end(); ++iElectron) {
            if (iElectron->pt() < 10) continue;

            TCElectron* eleCon = new ((*recoElectrons)[eleCount]) TCElectron;

            // Basic physics object info
            eleCon->SetPtEtaPhiM(iElectron->pt(), iElectron->eta(), iElectron->phi(), iElectron->mass());
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
            eleCon->SetSigmaIEtaIEta(iElectron->sigmaIetaIeta());
            eleCon->SetFBrem(iElectron->fbrem());
            eleCon->SetEOverP(iElectron->eSuperClusterOverP());
            eleCon->SetSCEta(iElectron->superCluster()->eta());
            eleCon->SetR9(iElectron->r9());

            eleCon->SetPtError(iElectron->gsfTrack()->ptError());
            eleCon->SetNormalizedChi2(iElectron->gsfTrack()->normalizedChi2());

            eleCon->SetNumberOfValidPixelHits(iElectron->gsfTrack()->hitPattern().numberOfValidPixelHits());
            eleCon->SetNumberOfValidTrackerHits(iElectron->gsfTrack()->hitPattern().numberOfValidTrackerHits());
            eleCon->SetNumberOfLostPixelHits(iElectron->gsfTrack()->hitPattern().numberOfLostPixelHits());
            eleCon->SetNumberOfLostTrackerHits(iElectron->gsfTrack()->hitPattern().numberOfLostTrackerHits());

            eleCon->SetIdMap("fabsEPDiff",fabs((1/iElectron->ecalEnergy()) - (1/iElectron->trackMomentumAtVtx().R()))); 

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
            eleCon->SetIsoMap("EffArea_R03", AEff03);
            eleCon->SetIsoMap("EffArea_R04", AEff04);

            // Conversion information
            bool convVeto = !(ConversionTools::hasMatchedConversion(*iElectron,hConversions,vertexBeamSpot.position()));
            eleCon->SetConversionVeto(convVeto);
            eleCon->SetConversionMissHits(iElectron->gsfTrack()->trackerExpectedHitsInner().numberOfHits());

            // Add electron MVA ID and ISO
            electronMVA(&(*iElectron), eleCon, iEvent, iSetup, thePfCollEleIso, rhoFactor);

            // Calculate electron energy regression
            InputTag  reducedEBRecHitCollection(string("reducedEcalRecHitsEB"));                                 
            InputTag  reducedEERecHitCollection(string("reducedEcalRecHitsEE"));                                 
            EcalClusterLazyTools lazyTools(iEvent, iSetup, reducedEBRecHitCollection, reducedEERecHitCollection);

            double eleEngReg = myEleReg->calculateRegressionEnergy(&(*iElectron), lazyTools, iSetup, rhoFactor, primaryVtcs->size(), false);
            double eleEngRegErr = myEleReg->calculateRegressionEnergyUncertainty(&(*iElectron), lazyTools, iSetup, rhoFactor, primaryVtcs->size(), false);


            eleCon->SetIdMap("EnergyRegression",eleEngReg);
            eleCon->SetIdMap("EnergyRegressionErr",eleEngRegErr);

            GsfElectron iElectronTmp = *iElectron;
            ElectronEnergyCalibrator myCalibrator("none",true,!isRealData,true,10,false);
            myCalibrator.correct(iElectronTmp, iElectronTmp.r9(), iEvent, iSetup, eleEngReg,eleEngRegErr);

            TLorentzVector tmpP4;
            tmpP4.SetPtEtaPhiE(iElectronTmp.pt(), iElectronTmp.eta(), iElectronTmp.phi(), iElectronTmp.energy());
            eleCon->SetRegressionMomCombP4(tmpP4);


            eleCount++;
        }
    }


    /////////////////
    // Get photons //
    /////////////////
    if (savePhotons_) {



        Handle<EcalRecHitCollection> Brechit;

        iEvent.getByLabel("reducedEcalRecHitsEB",Brechit);

        const EcalRecHitCollection* barrelRecHits= Brechit.product();

        Handle<vector<Photon> > photons;
        iEvent.getByLabel(photonTag_, photons);



        Handle<GsfElectronCollection> hElectrons;
        iEvent.getByLabel("gsfElectrons", hElectrons);


        for (vector<Photon>::const_iterator iPhoton = photons->begin(); iPhoton != photons->end() ; ++iPhoton) {


            TCPhoton* myPhoton = new ((*recoPhotons)[photonCount]) TCPhoton();

            //Crystal Info:
            std::vector< std::pair<DetId, float> >  PhotonHit_DetIds  = iPhoton->superCluster()->hitsAndFractions();
            std::vector<TCPhoton::CrystalInfo> crystalinfo_container;
            crystalinfo_container.clear();
            TCPhoton::CrystalInfo crystal = {};
            float timing_avg =0.0;
            int ncrys   = 0;
            int ncrysPhoton = 0;
            vector< std::pair<DetId, float> >::const_iterator detitr;

            for(detitr = PhotonHit_DetIds.begin(); detitr != PhotonHit_DetIds.end(); ++detitr)
            {

                if (((*detitr).first).det() == DetId::Ecal && ((*detitr).first).subdetId() == EcalBarrel) {
                    EcalRecHitCollection::const_iterator j= Brechit->find(((*detitr).first));
                    EcalRecHitCollection::const_iterator thishit;
                    if ( j!= Brechit->end())  thishit = j;
                    if ( j== Brechit->end()){
                        continue;
                    }

                    EBDetId detId  = (EBDetId)((*detitr).first);
                    crystal.rawId  = thishit->id().rawId();
                    crystal.energy = thishit->energy();
                    crystal.time   = thishit->time();
                    crystal.timeErr= thishit->timeError();
                    crystal.recoFlag = thishit->recoFlag();
                    crystal.ieta   = detId.ieta();
                    crystal.iphi   = detId.iphi();
                    if(crystal.energy > 0.1){
                        timing_avg  = timing_avg + crystal.time;
                        ncrys++;
                    }  
                }//end of if ((*detitr).det() == DetId::Ecal && (*detitr).subdetId() == EcalBarrel)
                crystalinfo_container.push_back(crystal);  
            }//End loop over detids
            std::sort(crystalinfo_container.begin(),crystalinfo_container.end(),EnergySortCriterium);


            //Without taking into account uncertainty, this time makes no sense.
            if (ncrys !=0) timing_avg = timing_avg/(float)ncrys;
            else timing_avg = -99.;
            ncrysPhoton = crystalinfo_container.size(); 
            int pho_timingavg_xtal      = timing_avg;

            myPhoton->SetNCrystals(crystalinfo_container.size());

            for (unsigned int y =0; y < crystalinfo_container.size() && y < 100;y++){ 

                myPhoton->SetCrystal(y,crystalinfo_container[y]);


            }//end of for (unsigned int y =0; y < crystalinfo_container.size();y++
            TCPhoton::CrystalInfo * savedCrystals = myPhoton->GetCrystalArray();
            /*
               for (int y = 0; y< myPhoton->GetNCrystals();y++){
               std::cout << "savedCrystals[y].time : " << savedCrystals[y].time << std::endl; 
               std::cout << "savedCrystals[y].timeErr : " << savedCrystals[y].timeErr << std::endl;
               std::cout << "savedCrystals[y].energy : " << savedCrystals[y].energy <<std::endl;
               std::cout << "savedCrystals[y].ieta: " << savedCrystals[y].ieta << std::endl;

               std::cout << "savedCrystals[y].rawId: " << savedCrystals[y].rawId <<std::endl;
               }
             */

            const BasicCluster& seedClus = *(iPhoton->superCluster()->seed());



            myPhoton->SetPxPyPzE(iPhoton->px(), iPhoton->py(), iPhoton->pz(), iPhoton->p());
            myPhoton->SetVtx(iPhoton->vx(), iPhoton->vy(), iPhoton->vz());

            // ID variables
            myPhoton->SetHadOverEm(iPhoton->hadTowOverEm());
            myPhoton->SetSigmaIEtaIEta(iPhoton->sigmaIetaIeta());
            myPhoton->SetR9(iPhoton->r9());
            myPhoton->SetTrackVeto(iPhoton->hasPixelSeed());

            myPhoton->SetSCEta(iPhoton->superCluster()->eta());
            myPhoton->SetSCPhi(iPhoton->superCluster()->phi());
            myPhoton->SetSCEnergy(iPhoton->superCluster()->energy());

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
        iEvent.getByLabel(InputTag("addPileupInfo"), PUInfo);
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

        if (saveGenParticles_) {
            Handle<GenParticleCollection> genParticleColl;
            iEvent.getByLabel("genParticles", genParticleColl);

            for (GenParticleCollection::const_iterator iGenPart = genParticleColl->begin(); iGenPart != genParticleColl->end(); ++iGenPart) {
                const GenParticle myParticle = GenParticle(*iGenPart);

                ////  Leptons and photons and b's, (oh my)
                if (
                        (
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
                        || abs(myParticle.pdgId()) == 35 
                        || abs(myParticle.pdgId()) == 36 
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

            Handle<GenJetCollection> GenJets;
            iEvent.getByLabel(genJetTag_, GenJets);

            for (GenJetCollection::const_iterator iJet = GenJets->begin(); iJet!= GenJets->end(); ++iJet) {
                GenJet myJet = GenJet(*iJet);      
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

    if (isRealData) {

        // Beam scraping filter
        myNoiseFilters.isScraping = isFilteredOutScraping(iEvent, iSetup, 10, 0.25);

        // ECAL filters
        Handle<bool> ecalTPFilterHandle;
        iEvent.getByLabel("EcalDeadCellTriggerPrimitiveFilter", ecalTPFilterHandle);
        if (ecalTPFilterHandle.isValid())  
            myNoiseFilters.isNoiseEcalTP = !(Bool_t)(*ecalTPFilterHandle);
        else 
            LogWarning("Filters")<<"Ecal TP NOT valid  ";

        Handle<bool> ecalBEFilterHandle;
        iEvent.getByLabel("EcalDeadCellBoundaryEnergyFilter", ecalBEFilterHandle);
        if (ecalBEFilterHandle.isValid())  
            myNoiseFilters.isNoiseEcalBE = !(Bool_t)(*ecalBEFilterHandle);
        else 
            LogWarning("Filters")<<"Ecal BE NOT valid  ";

        // HCAL
        Handle<bool> hcalNoiseFilterHandle;
        iEvent.getByLabel(hcalHBHEFilterTag_, hcalNoiseFilterHandle);
        if (hcalNoiseFilterHandle.isValid())  
            myNoiseFilters.isNoiseHcalHBHE = !(Bool_t)(*hcalNoiseFilterHandle);
        else 
            LogWarning("Filters")<<"hcal noise NOT valid  ";

        Handle<bool> hcalLaserFilterHandle;
        iEvent.getByLabel("hcalLaserEventFilter", hcalLaserFilterHandle);
        if (hcalLaserFilterHandle.isValid())  
            myNoiseFilters.isNoiseHcalLaser = !(Bool_t)(*hcalLaserFilterHandle);
        else 
            LogWarning("Filters")<<"hcal Laser NOT valid  ";

        Handle<bool> trackingFailureHandle;
        iEvent.getByLabel("trackingFailureFilter", trackingFailureHandle);
        if(trackingFailureHandle.isValid()) 
            myNoiseFilters.isNoiseTracking = !(Bool_t)(*trackingFailureHandle);
        else 
            LogWarning("Filters")<<"tracking Failure NOT valid  ";

        Handle<bool> eeBadScFilterHandle;
        iEvent.getByLabel("eeBadScFilter", eeBadScFilterHandle);
        if (eeBadScFilterHandle.isValid()) 
            myNoiseFilters.isNoiseEEBadSc = !(Bool_t)(*eeBadScFilterHandle);
        else 
            LogWarning("Filters")<<"eeBadSc NOT valid  ";


        Handle<bool> trkPOGFiltersHandle1;
        iEvent.getByLabel("manystripclus53X", trkPOGFiltersHandle1);
        if (trkPOGFiltersHandle1.isValid()) 
            myNoiseFilters.isNoisetrkPOG1 = !(Bool_t)(*trkPOGFiltersHandle1);
        else 
            LogWarning("Filters")<<"trkPOG1 NOT valid  ";

        Handle<bool> trkPOGFiltersHandle2;
        iEvent.getByLabel("toomanystripclus53X", trkPOGFiltersHandle2);
        if (trkPOGFiltersHandle2.isValid()) 
            myNoiseFilters.isNoisetrkPOG2 = !(Bool_t)(*trkPOGFiltersHandle2);
        else 
            LogWarning("Filters")<<"trkPOG2 NOT valid  ";

        Handle<bool> trkPOGFiltersHandle3;
        iEvent.getByLabel("logErrorTooManyClusters", trkPOGFiltersHandle3);
        if (trkPOGFiltersHandle3.isValid()) 
            myNoiseFilters.isNoisetrkPOG3 = !(Bool_t)(*trkPOGFiltersHandle3);
        else 
            LogWarning("Filters")<<"trkPOG3 NOT valid  ";


        Handle<BeamHaloSummary> TheBeamHaloSummary;
        iEvent.getByLabel("BeamHaloSummary",TheBeamHaloSummary);
        const BeamHaloSummary TheSummary = (*TheBeamHaloSummary.product() );
        if (!TheBeamHaloSummary.isValid()) 
            LogWarning("Filters")<<"The Summary (for CSC halo) NOT valid  ";

        myNoiseFilters.isCSCTightHalo = TheSummary.CSCTightHaloId();
        myNoiseFilters.isCSCLooseHalo = TheSummary.CSCLooseHaloId();


        //LogWarning("Filters")<<"\n csc1  "<< myNoiseFilters.isCSCTightHalo<<"  csc2  "<<myNoiseFilters.isCSCLooseHalo
        //  <<" isNoiseHcal HBHE "<<myNoiseFilters.isNoiseHcalHBHE<<"  laser "<<myNoiseFilters.isNoiseHcalLaser<<"\n"
        //  <<" ecal TP  "<<myNoiseFilters.isNoiseEcalTP<<"   ecal BE  "<<myNoiseFilters.isNoiseEcalBE;
    }

    ////////////////////////////  
    // get trigger information//
    ////////////////////////////

    Handle<TriggerResults> hltR;
    triggerResultsTag_ = InputTag(hlTriggerResults_,"",hltProcess_);
    iEvent.getByLabel(triggerResultsTag_,hltR);
    triggerEventTag_ = InputTag("hltTriggerSummaryAOD","",hltProcess_);
    Handle<trigger::TriggerEvent> hltE;                           
    iEvent.getByLabel(triggerEventTag_,hltE);                          

    const TriggerNames & triggerNames = iEvent.triggerNames(*hltR);
    hlNames = triggerNames.triggerNames();   

    triggerStatus   = ULong64_t(0x0);    

    for (int i=0; i < (int)hlNames.size(); ++i) {      
        if (!triggerDecision(hltR, i)) continue;	

        for (int j = 0; j < (int)triggerPaths_.size(); ++j){
            if (triggerPaths_[j] == "") continue;

            if (hlNames[i].compare(0, triggerPaths_[j].length(),triggerPaths_[j]) == 0) {
                //cout << hlNames[i] << " ?= " << triggerPaths_[j] << endl;
                triggerStatus |= ULong64_t(0x01) << j;

                /*
                   if (isRealData) {
                   pair<int, int> preScales;
                   preScales = hltConfig_.prescaleValues(iEvent, iSetup, hlNames[i]); 
                   hltPrescale[j] = preScales.first*preScales.second;
                   } else {
                 */
                hltPrescale[j] = 1;
                //}
            }
        }
    } 

    for(unsigned int t = 1; t<hlNames.size();t++){  
        analyzeTrigger(hltR, hltE, hlNames[t], &trigCount);       
    }                                               

    ++nEvents;

    /*if (eleCount == 0 || muCount == 0)*/  eventTree -> Fill(); // possibly specify a cut in configuration

    primaryVtx    -> Clear("C");
    recoJets      -> Clear("C");
    recoMuons     -> Clear("C");
    recoElectrons -> Clear("C");
    recoPhotons   -> Clear("C");
    triggerObjects-> Clear("C");
    genJets       -> Clear("C");
    genParticles  -> Clear("C");
}

// ------------ method called once each job just before starting event loop  ------------
void  ntupleProducer::beginJob()
{  
    eventTree      = fs->make<TTree>("eventTree","eventTree");
    jobTree        = fs->make<TTree>("jobTree", "jobTree");

    primaryVtx     = new TClonesArray("TCPrimaryVtx");
    recoJets       = new TClonesArray("TCJet");
    recoElectrons  = new TClonesArray("TCElectron");
    recoMuons      = new TClonesArray("TCMuon");
    recoPhotons    = new TClonesArray("TCPhoton");
    triggerObjects = new TClonesArray("TCTriggerObject");
    genJets        = new TClonesArray("TCGenJet");
    genParticles   = new TClonesArray("TCGenParticle");
    beamSpot       = new TVector3();
    recoMET        = 0;
    track_MET      = 0; 
    T0MET          = 0; 
    T2MET          = 0;
    mva_MET        = 0;

    eventTree->Branch("recoJets",&recoJets, 6400, 0);
    eventTree->Branch("recoElectrons",&recoElectrons, 6400, 0);
    eventTree->Branch("recoMuons",&recoMuons, 6400, 0);
    eventTree->Branch("recoPhotons",&recoPhotons, 6400, 0);
    eventTree->Branch("recoMET", &recoMET, 6400, 0);
    eventTree->Branch("mva_MET", &mva_MET, 6400, 0);
    eventTree->Branch("track_MET", &track_MET, 6400, 0);
    eventTree->Branch("T0MET", &T0MET, 6400, 0);
    eventTree->Branch("T2MET", &T2MET, 6400, 0);
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
    eventTree->Branch("rhoMuFactor",&rhoMuFactor, "rhoMuFactor/F");
    eventTree->Branch("triggerStatus",&triggerStatus, "triggerStatus/l");
    eventTree->Branch("hltPrescale",hltPrescale, "hltPrescale[64]/i");

    eventTree->Branch("NoiseFilters", &myNoiseFilters.isScraping, "isScraping/O:isNoiseHcalHBHE:isNoiseHcalLaser:isNoiseEcalTP:isNoiseEcalBE:isCSCTightHalo:isCSCLooseHalo:isNoiseTracking:isNoiseEEBadSc:isNoisetrkPOG1:isNoisetrkPOG2:isNoisetrkPOG3");

    jobTree->Branch("nEvents",&nEvents, "nEvents/i");
    jobTree->Branch("triggerNames", "vector<string>", &triggerPaths_);

    // Initialize HLT prescales //

    for (int i = 0; i < (int)(sizeof(hltPrescale)/sizeof(int)); ++i) hltPrescale[i] = 1;

    // Start counting number of events per job //
    nEvents = 0;

    // Photon Iso maker init
    phoIsolator.initializePhotonIsolation(kTRUE);
    phoIsolator.setConeSize(0.3);

    // Initialize Electron MVA nonsense
    eleIsolator.initializeElectronIsolation(kTRUE);
    eleIsolator.setConeSize(0.4);

    // Initialize Electron Regression
    myEleReg = new ElectronEnergyRegressionEvaluate();
    //myEleReg->initialize(mvaPath+"/src/data/eleEnergyRegWeights_V1.root",

    string mvaPath = getenv("CMSSW_BASE");
    mvaPath = mvaPath+"/src/EGamma/EGammaAnalysisTools/data:"+getenv("CMSSW_SEARCH_PATH");
    setenv("CMSSW_SEARCH_PATH",mvaPath.c_str(),1);

    myEleReg->initialize("eleEnergyRegWeights_V1.root", ElectronEnergyRegressionEvaluate::kNoTrkVar);

    if (verboseMVAs) cout<<"mvaPath: "<<mvaPath<<endl;
    if (verboseMVAs) cout<<"MVA electron regression shit probably has initialized"<<endl;
}

void ntupleProducer::beginRun(const Run& iRun, const EventSetup& iSetup)
{
    bool changed = true; 
    hltConfig_.init(iRun, iSetup, hltProcess_, changed);
    deliveredLumi = 0;
    recordedLumi  = 0;
}

void ntupleProducer::endLuminosityBlock(const LuminosityBlock& iLumi, const EventSetup& iSetup)
{
}


void ntupleProducer::endRun(const Run& iRun, const EventSetup& iSetup)
{
}


void ntupleProducer::endJob() 
{
    cout<<nEvents<<endl;
    jobTree->Fill();
}


bool ntupleProducer::triggerDecision(Handle<TriggerResults> &hltR, int iTrigger)
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


bool ntupleProducer::isFilteredOutScraping( const Event& iEvent, const EventSetup& iSetup, int numtrack, double thresh)
{

    bool  accepted = false;
    float fraction = 0;  
    // get GeneralTracks collection

    Handle<TrackCollection> tkRef;
    iEvent.getByLabel("generalTracks",tkRef);    
    const TrackCollection* tkColl = tkRef.product();

    int numhighpurity=0;
    TrackBase::TrackQuality _trackQuality = TrackBase::qualityByName("highPurity");

    if(tkColl->size()>(UInt_t)numtrack){ 
        TrackCollection::const_iterator itk = tkColl->begin();
        TrackCollection::const_iterator itk_e = tkColl->end();
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


bool ntupleProducer::associateJetToVertex(PFJet inJet, Handle<VertexCollection> vtxCollection, TCJet *outJet)
{
    if (fabs(inJet.eta()) > 2.5 || !inJet.getTrackRefs()) {
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
    vector<const Track*> jetTracks;
    float sumTrackX, sumTrackY, sumTrackZ, sumTrackPt;
    int   nJetTracks = 0;
    int   vCount = 0;

    sumTrackX = sumTrackY = sumTrackZ  = sumTrackPt = 0;

    const TrackRefVector &tracks = inJet.getTrackRefs(); 

    for (TrackRefVector::const_iterator iTrack = tracks.begin(); iTrack != tracks.end(); ++iTrack) {
        const Track &jetTrack = **iTrack;

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
            Vertex myVtx = Vertex(*iVtx); 
            if(!myVtx.isValid() || myVtx.isFake()) continue;
            associatedTrackSumPt.push_back(0);            
            associatedTrackCount.push_back(0);            

            for(Vertex::trackRef_iterator iTrackRef = myVtx.tracks_begin(); iTrackRef != myVtx.tracks_end(); ++iTrackRef){
                const RefToBase<Track> &myTrackRef = *iTrackRef; 

                if(myTrackRef.isAvailable()){
                    const Track &myVertexTrack = *myTrackRef.get();		

                    for(vector<const Track*>::const_iterator iTrack = jetTracks.begin(); iTrack != jetTracks.end(); ++iTrack){
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


void ntupleProducer::electronMVA(const GsfElectron* iElectron, TCElectron* eleCon, 
        const Event& iEvent, const EventSetup& iSetup, const PFCandidateCollection& PFCandidates, float Rho)
{
    if (verboseMVAs) cout<<"loading up electron MVA values"<<endl;
    //**********************************************************
    //ID variables
    //**********************************************************
    bool validKF= false; 
    TrackRef myTrackRef = iElectron->closestCtfTrackRef();
    validKF = (myTrackRef.isAvailable());
    validKF = (myTrackRef.isNonnull());  

    eleCon->SetIdMap("fbrem", (iElectron->fbrem() < -1) ? -1 : iElectron->fbrem());
    eleCon->SetIdMap("gsfChi2", (iElectron->gsfTrack()->normalizedChi2() > 200) ? 200 : iElectron->gsfTrack()->normalizedChi2());
    eleCon->SetIdMap("kfChi2", (validKF) ? myTrackRef->normalizedChi2() : 0);
    eleCon->SetIdMap("kfNLayers", (validKF) ? myTrackRef->hitPattern().trackerLayersWithMeasurement() : -1.);
    eleCon->SetIdMap("kfNLayersAll", (validKF) ? myTrackRef->numberOfValidHits() : -1.);
    eleCon->SetIdMap("dEta", (fabs(iElectron->deltaEtaSuperClusterTrackAtVtx()) > 0.06) ? 0.06 : fabs(iElectron->deltaEtaSuperClusterTrackAtVtx()));
    eleCon->SetIdMap("dPhi", iElectron->deltaPhiSuperClusterTrackAtVtx());
    eleCon->SetIdMap("dEtaAtCalo", iElectron->deltaEtaSeedClusterTrackAtCalo());
    eleCon->SetIdMap("SigmaIPhiIPhi", (isnan(iElectron->sigmaIphiIphi())) ? 0 : iElectron->sigmaIphiIphi() );
    eleCon->SetIdMap("SCEtaWidth", iElectron->superCluster()->etaWidth());
    eleCon->SetIdMap("SCPhiWidth", iElectron->superCluster()->phiWidth());
    eleCon->SetIdMap("EoP", (iElectron->eSuperClusterOverP() > 20) ? 20 : iElectron->eSuperClusterOverP());
    eleCon->SetIdMap("ome1x5oe5x5",(iElectron->e5x5()) !=0. ? 1.-(iElectron->e1x5()/iElectron->e5x5()) : -1.);
    if (eleCon->IdMap("ome1x5oe5x5") < -1) eleCon->SetIdMap("ome1x5oe5x5",-1);
    if (eleCon->IdMap("ome1x5oe5x5") > 2) eleCon->SetIdMap("ome1x5oe5x5",2);
    eleCon->SetIdMap("R9",(iElectron->r9() > 5) ? 5 : iElectron->r9());
    eleCon->SetIdMap("ooemoopV1",(1.0/iElectron->ecalEnergy()) - (1.0 / iElectron->p()));
    eleCon->SetIdMap("ooemoopV2",(1.0/iElectron->superCluster()->energy()) - (1.0 / iElectron->trackMomentumAtVtx().R()));
    eleCon->SetIdMap("eopOut",(iElectron->eEleClusterOverPout() > 20) ? 20 : iElectron->eEleClusterOverPout());
    eleCon->SetIdMap("preShowerORaw",iElectron->superCluster()->preshowerEnergy() / iElectron->superCluster()->rawEnergy());

    InputTag  vertexLabel(string("offlinePrimaryVertices"));
    Handle<VertexCollection> thePrimaryVertexColl;
    iEvent.getByLabel(vertexLabel,thePrimaryVertexColl);

    Vertex dummy;
    const Vertex *pv = &dummy;
    if (thePrimaryVertexColl->size() != 0) {
        pv = &*thePrimaryVertexColl->begin();
    } else { // create a dummy PV
        Vertex::Error e;
        e(0, 0) = 0.0015 * 0.0015;
        e(1, 1) = 0.0015 * 0.0015;
        e(2, 2) = 15. * 15.;
        Vertex::Point p(0, 0, 0);
        dummy = Vertex(p, e, 0, 0, 0);
    }

    //d0
    float fMVAVar_d0 = -9999.0;
    if (iElectron->gsfTrack().isNonnull()) {
        fMVAVar_d0 = (-1.0)*iElectron->gsfTrack()->dxy(pv->position()); 
    } else if (iElectron->closestCtfTrackRef().isNonnull()) {
        fMVAVar_d0 = (-1.0)*iElectron->closestCtfTrackRef()->dxy(pv->position()); 
    } else {
        fMVAVar_d0 = -9999.0;
    }

    eleCon->SetIdMap("d0",fMVAVar_d0);

    //default values for IP3D
    float fMVAVar_ip3d      = -999.0; 
    float fMVAVar_ip3dSig   = 0.0;

    ESHandle<TransientTrackBuilder> builder;
    iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", builder);
    TransientTrackBuilder thebuilder = *(builder.product());

    if (iElectron->gsfTrack().isNonnull()) {
        const double gsfsign   = ( (-iElectron->gsfTrack()->dxy(pv->position()))   >=0 ) ? 1. : -1.;

        const TransientTrack &tt = thebuilder.build(iElectron->gsfTrack()); 
        const std::pair<bool,Measurement1D> &ip3dpv =  IPTools::absoluteImpactParameter3D(tt,*pv);
        if (ip3dpv.first) {
            double ip3d = gsfsign*ip3dpv.second.value();
            double ip3derr = ip3dpv.second.error();  
            fMVAVar_ip3d = ip3d; 
            fMVAVar_ip3dSig = ip3d/ip3derr;
        }
    }

    eleCon->SetIdMap("ip3d",fMVAVar_ip3d);
    eleCon->SetIdMap("ip3dSig",fMVAVar_ip3dSig);

    //**********************************************************
    //Isolation variables
    //**********************************************************

    Double_t tmpChargedIso_DR0p0To0p1  = 0;
    Double_t tmpChargedIso_DR0p1To0p2  = 0;
    Double_t tmpChargedIso_DR0p2To0p3  = 0;
    Double_t tmpChargedIso_DR0p3To0p4  = 0;
    Double_t tmpChargedIso_DR0p4To0p5  = 0;
    Double_t tmpGammaIso_DR0p0To0p1  = 0;
    Double_t tmpGammaIso_DR0p1To0p2  = 0;
    Double_t tmpGammaIso_DR0p2To0p3  = 0;
    Double_t tmpGammaIso_DR0p3To0p4  = 0;
    Double_t tmpGammaIso_DR0p4To0p5  = 0;
    Double_t tmpNeutralHadronIso_DR0p0To0p1  = 0;
    Double_t tmpNeutralHadronIso_DR0p1To0p2  = 0;
    Double_t tmpNeutralHadronIso_DR0p2To0p3  = 0;
    Double_t tmpNeutralHadronIso_DR0p3To0p4  = 0;
    Double_t tmpNeutralHadronIso_DR0p4To0p5  = 0;

    //************************************************************
    //Note: Input collection is assumed to be PFNoPU collection
    //************************************************************
    for (PFCandidateCollection::const_iterator iP = PFCandidates.begin(); 
            iP != PFCandidates.end(); ++iP) {

        double dr = sqrt(pow(iP->eta() - iElectron->eta(),2) + pow(acos(cos(iP->phi() - iElectron->phi())),2));

        Bool_t passVeto = kTRUE;
        //Charged
        if(iP->trackRef().isNonnull()) {         

            //make sure charged pf candidates pass the PFNoPU condition (assumed)

            //************************************************************
            // Veto any PFmuon, or PFEle
            if (iP->particleId() == PFCandidate::e || iP->particleId() == PFCandidate::mu) passVeto = kFALSE;
            //************************************************************
            //************************************************************
            // Footprint Veto
            if (fabs(iElectron->superCluster()->eta()) > 1.479 && dr < 0.015) passVeto = kFALSE;
            if (iP->superClusterRef().isNonnull() && 
                    iP->superClusterRef() == iElectron->superCluster()) passVeto = kFALSE;
            if (iP->gsfTrackRef().isNonnull() && iElectron->gsfTrack().isNonnull() && 
                    iP->gsfTrackRef() == iElectron->gsfTrack()) passVeto = kFALSE;
            if (iP->trackRef().isNonnull() && iElectron->closestCtfTrackRef().isNonnull() && 
                    iP->trackRef() == iElectron->closestCtfTrackRef()) passVeto = kFALSE;
            //************************************************************
            if (passVeto) {
                if (dr < 0.1) tmpChargedIso_DR0p0To0p1 += iP->pt();
                if (dr >= 0.1 && dr < 0.2) tmpChargedIso_DR0p1To0p2 += iP->pt();
                if (dr >= 0.2 && dr < 0.3) tmpChargedIso_DR0p2To0p3 += iP->pt();
                if (dr >= 0.3 && dr < 0.4) tmpChargedIso_DR0p3To0p4 += iP->pt();
                if (dr >= 0.4 && dr < 0.5) tmpChargedIso_DR0p4To0p5 += iP->pt();
            } //pass veto    
        }
        //Gamma
        else if (iP->particleId() == PFCandidate::gamma) {
            //************************************************************
            // Footprint Veto
            if (fabs(iElectron->superCluster()->eta()) > 1.479 && dr < 0.08) passVeto = kFALSE;
            if (iP->superClusterRef() == iElectron->superCluster()) passVeto = kFALSE;
            //************************************************************  
            if (passVeto) {
                if (dr < 0.1) tmpGammaIso_DR0p0To0p1 += iP->pt();
                if (dr >= 0.1 && dr < 0.2) tmpGammaIso_DR0p1To0p2 += iP->pt();
                if (dr >= 0.2 && dr < 0.3) tmpGammaIso_DR0p2To0p3 += iP->pt();
                if (dr >= 0.3 && dr < 0.4) tmpGammaIso_DR0p3To0p4 += iP->pt();
                if (dr >= 0.4 && dr < 0.5) tmpGammaIso_DR0p4To0p5 += iP->pt();
            }
        }
        //NeutralHadron
        else {
            if (dr < 0.1) tmpNeutralHadronIso_DR0p0To0p1 += iP->pt();
            if (dr >= 0.1 && dr < 0.2) tmpNeutralHadronIso_DR0p1To0p2 += iP->pt();
            if (dr >= 0.2 && dr < 0.3) tmpNeutralHadronIso_DR0p2To0p3 += iP->pt();
            if (dr >= 0.3 && dr < 0.4) tmpNeutralHadronIso_DR0p3To0p4 += iP->pt();
            if (dr >= 0.4 && dr < 0.5) tmpNeutralHadronIso_DR0p4To0p5 += iP->pt();
        }
    } //loop over PF candidates

    double fMVAVar_ChargedIso_DR0p0To0p1,fMVAVar_ChargedIso_DR0p1To0p2,fMVAVar_ChargedIso_DR0p2To0p3,fMVAVar_ChargedIso_DR0p3To0p4,fMVAVar_ChargedIso_DR0p4To0p5,fMVAVar_GammaIso_DR0p0To0p1,fMVAVar_GammaIso_DR0p1To0p2,fMVAVar_GammaIso_DR0p2To0p3,
           fMVAVar_GammaIso_DR0p3To0p4,fMVAVar_GammaIso_DR0p4To0p5,fMVAVar_NeutralHadronIso_DR0p0To0p1,fMVAVar_NeutralHadronIso_DR0p1To0p2,fMVAVar_NeutralHadronIso_DR0p2To0p3,fMVAVar_NeutralHadronIso_DR0p3To0p4,fMVAVar_NeutralHadronIso_DR0p4To0p5;

    bool doPUCorrection = false;
    if (doPUCorrection) {
        fMVAVar_ChargedIso_DR0p0To0p1   = TMath::Min((tmpChargedIso_DR0p0To0p1)/iElectron->pt(), 2.5);
        fMVAVar_ChargedIso_DR0p1To0p2   = TMath::Min((tmpChargedIso_DR0p1To0p2)/iElectron->pt(), 2.5);
        fMVAVar_ChargedIso_DR0p2To0p3   = TMath::Min((tmpChargedIso_DR0p2To0p3)/iElectron->pt(), 2.5);
        fMVAVar_ChargedIso_DR0p3To0p4   = TMath::Min((tmpChargedIso_DR0p3To0p4)/iElectron->pt(), 2.5);
        fMVAVar_ChargedIso_DR0p4To0p5   = TMath::Min((tmpChargedIso_DR0p4To0p5)/iElectron->pt(), 2.5); 
        fMVAVar_GammaIso_DR0p0To0p1     = TMath::Max(TMath::Min((tmpGammaIso_DR0p0To0p1 - Rho*ElectronEffectiveArea::GetElectronEffectiveArea(ElectronEffectiveArea::kEleGammaIsoDR0p0To0p1, iElectron->superCluster()->eta(), ElectronEffectiveArea::kEleEAData2012))/iElectron->pt(), 2.5), 0.0);
        fMVAVar_GammaIso_DR0p1To0p2     = TMath::Max(TMath::Min((tmpGammaIso_DR0p1To0p2 - Rho*ElectronEffectiveArea::GetElectronEffectiveArea(ElectronEffectiveArea::kEleGammaIsoDR0p1To0p2, iElectron->superCluster()->eta(), ElectronEffectiveArea::kEleEAData2012))/iElectron->pt(), 2.5), 0.0);
        fMVAVar_GammaIso_DR0p2To0p3 = TMath::Max(TMath::Min((tmpGammaIso_DR0p2To0p3 - Rho*ElectronEffectiveArea::GetElectronEffectiveArea(ElectronEffectiveArea::kEleGammaIsoDR0p2To0p3, iElectron->superCluster()->eta(), ElectronEffectiveArea::kEleEAData2012))/iElectron->pt(), 2.5), 0.0);
        fMVAVar_GammaIso_DR0p3To0p4 = TMath::Max(TMath::Min((tmpGammaIso_DR0p3To0p4 - Rho*ElectronEffectiveArea::GetElectronEffectiveArea(ElectronEffectiveArea::kEleGammaIsoDR0p3To0p4, iElectron->superCluster()->eta(), ElectronEffectiveArea::kEleEAData2012))/iElectron->pt(), 2.5), 0.0);
        fMVAVar_GammaIso_DR0p4To0p5 = TMath::Max(TMath::Min((tmpGammaIso_DR0p4To0p5 - Rho*ElectronEffectiveArea::GetElectronEffectiveArea(ElectronEffectiveArea::kEleGammaIsoDR0p4To0p5, iElectron->superCluster()->eta(), ElectronEffectiveArea::kEleEAData2012))/iElectron->pt(), 2.5), 0.0);
        fMVAVar_NeutralHadronIso_DR0p0To0p1 = TMath::Max(TMath::Min((tmpNeutralHadronIso_DR0p0To0p1 - Rho*ElectronEffectiveArea::GetElectronEffectiveArea(ElectronEffectiveArea::kEleNeutralHadronIsoDR0p0To0p1, iElectron->superCluster()->eta(), ElectronEffectiveArea::kEleEAData2012))/iElectron->pt(), 2.5), 0.0);
        fMVAVar_NeutralHadronIso_DR0p1To0p2 = TMath::Max(TMath::Min((tmpNeutralHadronIso_DR0p1To0p2 - Rho*ElectronEffectiveArea::GetElectronEffectiveArea(ElectronEffectiveArea::kEleNeutralHadronIsoDR0p1To0p2, iElectron->superCluster()->eta(), ElectronEffectiveArea::kEleEAData2012))/iElectron->pt(), 2.5), 0.0);
        fMVAVar_NeutralHadronIso_DR0p2To0p3 = TMath::Max(TMath::Min((tmpNeutralHadronIso_DR0p2To0p3 - Rho*ElectronEffectiveArea::GetElectronEffectiveArea(ElectronEffectiveArea::kEleNeutralHadronIsoDR0p2To0p3, iElectron->superCluster()->eta(), ElectronEffectiveArea::kEleEAData2012))/iElectron->pt(), 2.5), 0.0);
        fMVAVar_NeutralHadronIso_DR0p3To0p4 = TMath::Max(TMath::Min((tmpNeutralHadronIso_DR0p3To0p4 - Rho*ElectronEffectiveArea::GetElectronEffectiveArea(ElectronEffectiveArea::kEleNeutralHadronIsoDR0p3To0p4, iElectron->superCluster()->eta(), ElectronEffectiveArea::kEleEAData2012))/iElectron->pt(), 2.5), 0.0);
        fMVAVar_NeutralHadronIso_DR0p4To0p5 = TMath::Max(TMath::Min((tmpNeutralHadronIso_DR0p4To0p5 - Rho*ElectronEffectiveArea::GetElectronEffectiveArea(ElectronEffectiveArea::kEleNeutralHadronIsoDR0p4To0p5, iElectron->superCluster()->eta(), ElectronEffectiveArea::kEleEAData2012))/iElectron->pt(), 2.5), 0.0);
    } else {
        fMVAVar_ChargedIso_DR0p0To0p1   = TMath::Min((tmpChargedIso_DR0p0To0p1)/iElectron->pt(), 2.5);
        fMVAVar_ChargedIso_DR0p1To0p2   = TMath::Min((tmpChargedIso_DR0p1To0p2)/iElectron->pt(), 2.5) / 0.03;
        fMVAVar_ChargedIso_DR0p2To0p3 = TMath::Min((tmpChargedIso_DR0p2To0p3)/iElectron->pt(), 2.5) / 0.05;
        fMVAVar_ChargedIso_DR0p3To0p4 = TMath::Min((tmpChargedIso_DR0p3To0p4)/iElectron->pt(), 2.5) / 0.07;
        fMVAVar_ChargedIso_DR0p4To0p5 = TMath::Min((tmpChargedIso_DR0p4To0p5)/iElectron->pt(), 2.5) / 0.09; 
        fMVAVar_GammaIso_DR0p0To0p1 = TMath::Max(TMath::Min((tmpGammaIso_DR0p0To0p1)/iElectron->pt(), 2.5), 0.0);
        fMVAVar_GammaIso_DR0p1To0p2 = TMath::Max(TMath::Min((tmpGammaIso_DR0p1To0p2)/iElectron->pt(), 2.5), 0.0) / 0.03;
        fMVAVar_GammaIso_DR0p2To0p3 = TMath::Max(TMath::Min((tmpGammaIso_DR0p2To0p3)/iElectron->pt(), 2.5), 0.0) / 0.05;
        fMVAVar_GammaIso_DR0p3To0p4 = TMath::Max(TMath::Min((tmpGammaIso_DR0p3To0p4)/iElectron->pt(), 2.5), 0.0) / 0.07;
        fMVAVar_GammaIso_DR0p4To0p5 = TMath::Max(TMath::Min((tmpGammaIso_DR0p4To0p5)/iElectron->pt(), 2.5), 0.0) / 0.09;
        fMVAVar_NeutralHadronIso_DR0p0To0p1 = TMath::Max(TMath::Min((tmpNeutralHadronIso_DR0p0To0p1)/iElectron->pt(), 2.5), 0.0);
        fMVAVar_NeutralHadronIso_DR0p1To0p2 = TMath::Max(TMath::Min((tmpNeutralHadronIso_DR0p1To0p2)/iElectron->pt(), 2.5), 0.0) / 0.03;
        fMVAVar_NeutralHadronIso_DR0p2To0p3 = TMath::Max(TMath::Min((tmpNeutralHadronIso_DR0p2To0p3)/iElectron->pt(), 2.5), 0.0) / 0.05;
        fMVAVar_NeutralHadronIso_DR0p3To0p4 = TMath::Max(TMath::Min((tmpNeutralHadronIso_DR0p3To0p4)/iElectron->pt(), 2.5), 0.0) / 0.07;
        fMVAVar_NeutralHadronIso_DR0p4To0p5 = TMath::Max(TMath::Min((tmpNeutralHadronIso_DR0p4To0p5)/iElectron->pt(), 2.5), 0.0) / 0.09;
    }

    eleCon->SetIsoMap("ChargedIso_DR0p0To0p1",fMVAVar_ChargedIso_DR0p0To0p1);
    eleCon->SetIsoMap("ChargedIso_DR0p1To0p2",fMVAVar_ChargedIso_DR0p1To0p2);
    eleCon->SetIsoMap("ChargedIso_DR0p2To0p3",fMVAVar_ChargedIso_DR0p2To0p3);
    eleCon->SetIsoMap("ChargedIso_DR0p3To0p4",fMVAVar_ChargedIso_DR0p3To0p4);
    eleCon->SetIsoMap("ChargedIso_DR0p4To0p5",fMVAVar_ChargedIso_DR0p4To0p5);
    eleCon->SetIsoMap("GammaIso_DR0p0To0p1",fMVAVar_GammaIso_DR0p0To0p1);
    eleCon->SetIsoMap("GammaIso_DR0p1To0p2",fMVAVar_GammaIso_DR0p1To0p2);
    eleCon->SetIsoMap("GammaIso_DR0p2To0p3",fMVAVar_GammaIso_DR0p2To0p3);
    eleCon->SetIsoMap("GammaIso_DR0p3To0p4",fMVAVar_GammaIso_DR0p3To0p4);
    eleCon->SetIsoMap("GammaIso_DR0p4To0p5",fMVAVar_GammaIso_DR0p4To0p5);
    eleCon->SetIsoMap("NeutralHadronIso_DR0p0To0p1",fMVAVar_NeutralHadronIso_DR0p0To0p1);
    eleCon->SetIsoMap("NeutralHadronIso_DR0p1To0p2",fMVAVar_NeutralHadronIso_DR0p1To0p2);
    eleCon->SetIsoMap("NeutralHadronIso_DR0p2To0p3",fMVAVar_NeutralHadronIso_DR0p2To0p3);
    eleCon->SetIsoMap("NeutralHadronIso_DR0p3To0p4",fMVAVar_NeutralHadronIso_DR0p3To0p4);
    eleCon->SetIsoMap("NeutralHadronIso_DR0p4To0p5",fMVAVar_NeutralHadronIso_DR0p4To0p5);

    bool preSelPassV2 = false;

    double electronTrackZ = 0;
    if (iElectron->gsfTrack().isNonnull()) {
        electronTrackZ = iElectron->gsfTrack()->dz(pv->position());
    } else if (iElectron->closestCtfTrackRef().isNonnull()) {
        electronTrackZ = iElectron->closestCtfTrackRef()->dz(pv->position());
    }

    if (fabs(iElectron->superCluster()->eta()) < 1.479) {
        if ( (
                    iElectron->sigmaIetaIeta()< 0.01
                    && fabs(iElectron->deltaEtaSuperClusterTrackAtVtx()) < 0.007
                    && fabs(iElectron->deltaPhiSuperClusterTrackAtVtx()) < 0.15
                    && iElectron->hadronicOverEm() < 0.12
                    && iElectron->gsfTrack()->trackerExpectedHitsInner().numberOfLostHits() <= 1
                    && fabs(electronTrackZ) < 0.1
                    && ( iElectron->dr03TkSumPt()) / iElectron->pt() < 0.2
                    && ( fmax(iElectron->dr03EcalRecHitSumEt() - 1.0,0.0) ) /iElectron->pt() < 0.20
                    && (iElectron->dr03HcalTowerSumEt()) / iElectron->pt() < 0.20
             )
           ) {
            preSelPassV2= true;
        }
    } else { //endcap
        if ( (  
                    iElectron->sigmaIetaIeta()< 0.03
                    && fabs(iElectron->deltaEtaSuperClusterTrackAtVtx()) < 0.009
                    && fabs(iElectron->deltaPhiSuperClusterTrackAtVtx()) < 0.10
                    && iElectron->hadronicOverEm() < 0.10
                    && iElectron->gsfTrack()->trackerExpectedHitsInner().numberOfLostHits() <= 1
                    && fabs(electronTrackZ) < 0.1
                    && (iElectron->dr03TkSumPt() ) / iElectron->pt() < 0.2
                    && (iElectron->dr03EcalRecHitSumEt() ) / iElectron->pt() < 0.20
                    && (iElectron->dr03HcalTowerSumEt()) / iElectron->pt() < 0.20

             )
           ) {
            preSelPassV2 = true;
        }
    }
    eleCon->SetIdMap("preSelPassV2",preSelPassV2);

    bool preSelPassV1 = false;

    if (fabs(iElectron->superCluster()->eta()) < 1.479) {
        if ( (
                    iElectron->sigmaIetaIeta()< 0.014
                    && iElectron->hadronicOverEm() < 0.15
                    && iElectron->gsfTrack()->trackerExpectedHitsInner().numberOfLostHits() == 0
                    && ( iElectron->dr03TkSumPt()) / iElectron->pt() < 0.2
                    && ( iElectron->dr03EcalRecHitSumEt()) /iElectron->pt() < 0.2
                    && (iElectron->dr03HcalTowerSumEt()) / iElectron->pt() < 0.2
             )
           ) {
            preSelPassV1 = true;
        }
    } else { //endcap
        if ( (  
                    iElectron->sigmaIetaIeta()< 0.035
                    && iElectron->hadronicOverEm() < 0.10
                    && iElectron->gsfTrack()->trackerExpectedHitsInner().numberOfLostHits() == 0
                    && (iElectron->dr03TkSumPt() ) / iElectron->pt() < 0.2
                    && (iElectron->dr03EcalRecHitSumEt() ) / iElectron->pt() < 0.20
                    && (iElectron->dr03HcalTowerSumEt()) / iElectron->pt() < 0.20

             )
           ) {
            preSelPassV1 = true;
        }
    }

    eleCon->SetIdMap("preSelPassV1",preSelPassV1);

    bool preSelPassV3 = false;
    if (fabs(iElectron->superCluster()->eta()) < 1.479) {
        if ( (
                    iElectron->sigmaIetaIeta()< 0.014
                    && iElectron->hadronicOverEm() < 0.15
                    && iElectron->gsfTrack()->trackerExpectedHitsInner().numberOfLostHits() == 0
                    && ( iElectron->dr03TkSumPt()) / iElectron->pt() < 0.2
                    && ( iElectron->dr03EcalRecHitSumEt()) /iElectron->pt() < 0.2
                    && (iElectron->dr03HcalTowerSumEt()) / iElectron->pt() < 0.2
             )
           ) {
            preSelPassV3 = true;
        }
    } else { //endcap
        if ( (  
                    iElectron->sigmaIetaIeta()< 0.035
                    && iElectron->hadronicOverEm() < 0.10
                    && iElectron->gsfTrack()->trackerExpectedHitsInner().numberOfLostHits() == 0
                    && (iElectron->dr03TkSumPt() ) / iElectron->pt() < 0.2
                    && (iElectron->dr03EcalRecHitSumEt()-1.0 ) / iElectron->pt() < 0.20
                    && (iElectron->dr03HcalTowerSumEt()) / iElectron->pt() < 0.20

             )
           ) {
            preSelPassV3 = true;
        }
    }
    eleCon->SetIdMap("preSelPassV3",preSelPassV3);

    return;
}

void ntupleProducer::analyzeTrigger(Handle<TriggerResults> &hltR, Handle<trigger::TriggerEvent> &hltE, const std::string& triggerName, int* trigCount) {

    using namespace trigger;

    const unsigned int n(hltConfig_.size());
    const unsigned int triggerIndex(hltConfig_.triggerIndex(triggerName));

    TLorentzVector triggerLepton;

    // abort on invalid trigger name

    bool goodTrigger = false;
    for (unsigned int i =0; i< triggerPaths_.size(); i++){
        if (triggerName.find(triggerPaths_[i]) != string::npos){
            goodTrigger = true;
            break;
        }
    }

    if(!goodTrigger) return;

    //if(verboseTrigs){
    //  std::cout<<" n = "<<n<<" triggerIndex = "<<triggerIndex<<" size = "<<hltConfig_.size()<<std::endl;
    //  std::cout<<" Analyze triggerName : "<<triggerName<<std::endl;
    //}
    //std::cout<<" Analyze triggerName : "<<triggerName<<std::endl;

    if (triggerIndex>=n) {
        if(verboseTrigs){
            cout << "DimuonAna::analyzeTrigger: path "
                << triggerName << " - not found!" << endl;
        }
        return;
    }

    // modules on this trigger path
    // const unsigned int moduleIndex(hltR->index(triggerIndex));
    const unsigned int moduleIndex(hltR->index(triggerIndex));
    const unsigned int m(hltConfig_.size(triggerIndex));
    const vector<string>& moduleLabels(hltConfig_.moduleLabels(triggerIndex));
    if (moduleIndex != m-1) return;
    if(verboseTrigs){
        cout << "DimuonAna::analyzeTrigger: path "
            << triggerName << " [" << triggerIndex << "]" << endl;

        std::cout<<"  n = "<< n<<" triggerIndex = "<<triggerIndex<<" m = "<<m<<std::endl;
        std::cout<<" moduleLabels = "<<moduleLabels.size()<<" moduleIndex = "<<moduleIndex<<std::endl;

        // Results from TriggerResults product
        cout << " Trigger path status:"
            << " WasRun=" << hltR->wasrun(triggerIndex)
            << " Accept=" << hltR->accept(triggerIndex)
            << " Error =" << hltR->error(triggerIndex)
            << endl;
        cout << " Last active module - label/type: "
            << moduleLabels[moduleIndex] << "/" << hltConfig_.moduleType(moduleLabels[moduleIndex])
            << " [" << moduleIndex << " out of 0-" << (m-1) << " on this path]"
            << endl;
    }
    assert (moduleIndex<m);

    // Results from TriggerEvent product - Attention: must look only for
    // modules actually run in this path for this event!
    std::vector < GlobalVector > passMomenta; 
    for (unsigned int j=0; j<=moduleIndex; ++j) {
        const string& moduleLabel(moduleLabels[j]);
        const string  moduleType(hltConfig_.moduleType(moduleLabel));

        // check whether the module is packed up in TriggerEvent product
        //cout<<hltE->filterIndex(InputTag(moduleLabel,"",hltProcess_))<<endl;

        const unsigned int filterIndex(hltE->filterIndex(InputTag(moduleLabel,"",hltProcess_)));

        //  if ( (moduleLabel.find("Calo") == string::npos) )continue;
        //  if ( (moduleLabel.find("hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoDZ") == string::npos)
        //      && (moduleLabel.find("hltEle17CaloId") == string::npos)
        //      && (moduleLabel.find("hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoDoubleFilter") == string::npos) ) continue;

        if(verboseTrigs){
            std::cout<<" j = "<<j<<" modLabel/moduleType = "<<moduleLabel<<"/"<<moduleType<<" filterIndex = "<<filterIndex<<" sizeF = "<<hltE->sizeFilters()<<std::endl;
        }
        if (filterIndex<hltE->sizeFilters()) {
            if(verboseTrigs){
                cout << " 'L3' (or 'L1', 'L2') filter in slot " << j << " - label/type " << moduleLabel << "/" << moduleType << endl;
            }
            const Vids& VIDS (hltE->filterIds(filterIndex));
            const Keys& KEYS(hltE->filterKeys(filterIndex));
            const size_type nI(VIDS.size());
            const size_type nK(KEYS.size());
            assert(nI==nK);
            const size_type n(max(nI,nK));
            if(verboseTrigs){
                cout << "   " << n  << " accepted 'L3' (or 'L1', 'L2') objects found: " << endl;
            }
            const TriggerObjectCollection& TOC(hltE->getObjects());
            for (size_type i=0; i!=n; ++i) {
                if(0==i){
                    passMomenta.clear();
                }
                const TriggerObject& TO(TOC[KEYS[i]]);
                GlobalVector momentumT0(TO.px(),TO.py(),TO.pz());
                if (TO.pt() < 10) continue;
                TCTriggerObject* trigObj = new ((*triggerObjects)[*trigCount]) TCTriggerObject;

                //std::cout<<" i_KEY = "<<i<<" id = "<<TO.id()<<" typ = "<<moduleType<<std::endl;
                //if("HLTLevel1GTSeed"==moduleType){
                //allMuL1TriggerVectors.push_back(momentumT0);
                ////std::cout<<" L1 object found"<<std::endl;
                //}

                if(verboseTrigs){
                    std::cout<<" i = "<<i<<" moduleLabel/moduleType : "<<moduleLabel<<"/"<<moduleType<<std::endl;
                    cout << "   " << i << " " << VIDS[i] << "/" << KEYS[i] << ": "
                        << TO.id() << " " << TO.pt() << " " << TO.eta() << " " << TO.phi() << " " << TO.mass()
                        << endl;
                }

                trigObj->SetPtEtaPhiE(TO.pt(),TO.eta(),TO.phi(),TO.energy());
                trigObj->SetHLTName(triggerName);
                trigObj->SetModuleName(moduleLabel);
                trigObj->SetId(TO.id());

                (*trigCount)++;

            }
            //
        }
    }
    //cout<<endl;
    return;
}

float ntupleProducer::MatchBTagsToJets(const JetTagCollection bTags,const PFJet jet)
{
    float discrValue = -999;
    for (JetTagCollection::const_iterator iTag = bTags.begin(); iTag != bTags.end(); iTag++) {
        if (sqrt(pow(iTag->first->eta() - jet.eta(), 2) + pow(deltaPhi(iTag->first->phi(),jet.phi()), 2)) == 0.) {
            discrValue = iTag->second;
            break;
        }
    }

    return discrValue;
}

//define this as a plug-in
DEFINE_FWK_MODULE(ntupleProducer);

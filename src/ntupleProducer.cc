#include "../interface/ntupleProducer.h"

ntupleProducer::ntupleProducer(const edm::ParameterSet& iConfig):
    eventTree(0)
{
    jetTag_           = iConfig.getUntrackedParameter<edm::InputTag>("JetTag");
    jecTag_           = iConfig.getParameter<std::string>("JecTag");
    metTag_           = iConfig.getUntrackedParameter<edm::InputTag>("METTag");
    muonTag_          = iConfig.getUntrackedParameter<edm::InputTag>("MuonTag");
    electronTag_      = iConfig.getUntrackedParameter<edm::InputTag>("ElectronTag");
    photonTag_        = iConfig.getUntrackedParameter<edm::InputTag>("PhotonTag");
    genJetTag_        = iConfig.getUntrackedParameter<edm::InputTag>("GenJetTag");
    primaryVtxTag_    = iConfig.getUntrackedParameter<edm::InputTag>("PrimaryVtxTag");

    rhoCorrTag_       = iConfig.getUntrackedParameter<edm::InputTag>("rhoCorrTag");
    rho25CorrTag_     = iConfig.getUntrackedParameter<edm::InputTag>("rho25CorrTag");
    rhoMuCorrTag_     = iConfig.getUntrackedParameter<edm::InputTag>("rhoMuCorrTag");
    hlTriggerResults_ = iConfig.getUntrackedParameter<string>("HLTriggerResults","TriggerResults");
    hltProcess_       = iConfig.getUntrackedParameter<string>("hltName");
    triggerPaths_     = iConfig.getUntrackedParameter<vector<string> >("triggers");

    partFlowTag_      = iConfig.getUntrackedParameter<edm::InputTag>("partFlowTag");

    saveEleCrystals_  = iConfig.getUntrackedParameter<bool>("saveEleCrystals");
    savePhoCrystals_  = iConfig.getUntrackedParameter<bool>("savePhoCrystals");

    saveTriggerObj_   = iConfig.getUntrackedParameter<bool>("saveTriggerObj");

    verboseTrigs       = iConfig.getUntrackedParameter<bool>("verboseTrigs");
    verboseMVAs        = iConfig.getUntrackedParameter<bool>("verboseMVAs");

    ecalTPFilterTag_    = iConfig.getUntrackedParameter<edm::InputTag>("ecalTPFilterTag");
    ecalBEFilterTag_    = iConfig.getUntrackedParameter<edm::InputTag>("ecalBEFilterTag");
    hcalHBHEFilterTag_  = iConfig.getUntrackedParameter<edm::InputTag>("hcalHBHEFilterTag");
    hcalLaserFilterTag_ = iConfig.getUntrackedParameter<edm::InputTag>("hcalLaserFilterTag");
    trackingFailureTag_ = iConfig.getUntrackedParameter<edm::InputTag>("trackingFailureTag");
    eeBadScFilterTag_   = iConfig.getUntrackedParameter<edm::InputTag>("eeBadScFilterTag");
    trkPOGFiltersTag1_  = iConfig.getUntrackedParameter<edm::InputTag>("trkPOGFiltersTag1");
    trkPOGFiltersTag2_  = iConfig.getUntrackedParameter<edm::InputTag>("trkPOGFiltersTag2");
    trkPOGFiltersTag3_  = iConfig.getUntrackedParameter<edm::InputTag>("trkPOGFiltersTag3");
    photonIsoCalcTag_   = iConfig.getParameter<edm::ParameterSet>("photonIsoCalcTag");
    jetPUIdAlgo_        = iConfig.getParameter<edm::ParameterSet>("jetPUIdAlgo");

    SCFPRemovalConePho_    = iConfig.getUntrackedParameter<double>("isolation_cone_size_forSCremoval_pho");
    SCFPRemovalConeEl_     = iConfig.getUntrackedParameter<double>("isolation_cone_size_forSCremoval_el");

    ebReducedRecHitCollection_ = iConfig.getParameter<edm::InputTag>("ebReducedRecHitCollection");
    eeReducedRecHitCollection_ = iConfig.getParameter<edm::InputTag>("eeReducedRecHitCollection");
    esReducedRecHitCollection_ = iConfig.getParameter<edm::InputTag>("esReducedRecHitCollection");

    // CiCPhoton class
    cicPhotonId_.reset(new CiCPhotonID(iConfig));
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
    // this is for CVS check test
    eventNumber  = iEvent.id().event();
    runNumber    = iEvent.id().run();
    lumiSection  = (unsigned int)iEvent.getLuminosityBlock().luminosityBlock();
    bunchCross   = (unsigned int)iEvent.bunchCrossing();
    isRealData   = iEvent.isRealData();

    edm::Handle<reco::BeamSpot> beamSpotHandle;
    iEvent.getByLabel("offlineBeamSpot", beamSpotHandle);
    reco::BeamSpot vertexBeamSpot = *beamSpotHandle;

    beamSpot->SetXYZ(vertexBeamSpot.x0(), vertexBeamSpot.y0(), vertexBeamSpot.z0());

    int vtxCount, jetCount, metCount, muCount, pfMuCount, eleCount, photonCount, pfPhotonCount, genCount, genPartCount;
    vtxCount = jetCount = metCount = muCount = pfMuCount = eleCount = photonCount = pfPhotonCount = genCount = genPartCount = 0;


    /////////////////////////////////////
    // Get PF candidates for later use //
    /////////////////////////////////////


    Handle<PFCandidateCollection> pfCands;
    iEvent.getByLabel(partFlowTag_,pfCands);
    const  PFCandidateCollection thePfColl = *(pfCands.product());

    //get a lazyTool

    lazyTool.reset(new EcalClusterLazyTools(iEvent, iSetup, ebReducedRecHitCollection_, eeReducedRecHitCollection_));

    //configure CiCPhoton

    Handle<reco::VertexCollection> recVtxsBS_;
    iEvent.getByLabel(edm::InputTag("offlinePrimaryVerticesWithBS",""), recVtxsBS_);

    Handle<reco::TrackCollection> tracksHandle_;
    iEvent.getByLabel("generalTracks",tracksHandle_);

    Handle<reco::GsfElectronCollection > gsfElectronHandle_;
    iEvent.getByLabel(electronTag_, gsfElectronHandle_);

    Handle<double> rhoCorr;
    iEvent.getByLabel(rhoCorrTag_, rhoCorr);
    rhoFactor = (float)(*rhoCorr);


    cicPhotonId_->configure(recVtxsBS_, tracksHandle_, gsfElectronHandle_, pfCands, rhoFactor);

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

    VertexRef myVtxRef(primaryVtcs, 0); // main vertex #0

    ////////////////////////////
    // get trigger information//
    ////////////////////////////

    edm::Handle<TriggerResults> hltResults;
    triggerResultsTag_ = InputTag(hlTriggerResults_,"",hltProcess_);
    iEvent.getByLabel(triggerResultsTag_,hltResults);

    edm::Handle<trigger::TriggerEvent> hltEvent;
    triggerEventTag_ = InputTag("hltTriggerSummaryAOD","",hltProcess_);
    iEvent.getByLabel(triggerEventTag_,hltEvent);

    const TriggerNames & triggerNames = iEvent.triggerNames(*hltResults);
    hlNames = triggerNames.triggerNames();

    triggerStatus   = ULong64_t(0x0);

    int trigCount = 0;
    for (int i=0; i < (int)hlNames.size(); ++i) {
        if (!triggerDecision(hltResults, i)) continue;

        for (int j = 0; j < (int)triggerPaths_.size(); ++j){
            if (triggerPaths_[j] == "") continue;

            if (hlNames[i].compare(0, triggerPaths_[j].length(),triggerPaths_[j]) == 0) {
                //cout <<eventNumber<<"\n \t"<< hlNames[i] << " ?= " << triggerPaths_[j] <<"   is PASSED"<< endl;
                triggerStatus |= ULong64_t(0x01) << j;
                hltPrescale[j] = 1;

                if(saveTriggerObj_) analyzeTrigger(hltResults, hltEvent, hlNames[i], &trigCount);
                /* if (isRealData) {
                   pair<int, int> preScales;
                   preScales = hltConfig_.prescaleValues(iEvent, iSetup, hlNames[i]);
                   hltPrescale[j] = preScales.first*preScales.second;
                   } */
            }
        }
    }

    //for (unsigned i = 0; i < triggerObjects.size(); ++i) 
    //    if (fabs(triggerObjects[i].GetId()) == 13)
    //        cout << triggerObjects[i].Pt() << ", " << triggerObjects[i].Eta() << endl;


    ///////////////////////
    //get jet information//
    ///////////////////////

    Handle<double> rho25Corr;
    iEvent.getByLabel(rho25CorrTag_, rho25Corr);
    rho25Factor = (float)(*rho25Corr);

    Handle<double> rhoMuCorr;
    iEvent.getByLabel(rhoMuCorrTag_, rhoMuCorr);
    rhoMuFactor = (float)(*rhoMuCorr);

    //cout<<" RHOS. In eta 4.4 = "<<rhoFactor<<"   in eta25 "<<rho25Factor<<"  MUs: "<<rhoMuFactor<<endl;


    edm::Handle<reco::JetTagCollection> bTagCollectionTCHE;
    iEvent.getByLabel("trackCountingHighEffBJetTags", bTagCollectionTCHE);
    const reco::JetTagCollection & bTagsTCHE = *(bTagCollectionTCHE.product());

    edm::Handle<reco::JetTagCollection> bTagCollectionTCHP;
    iEvent.getByLabel("trackCountingHighPurBJetTags", bTagCollectionTCHP);
    const reco::JetTagCollection & bTagsTCHP = *(bTagCollectionTCHP.product());

    edm::Handle<reco::JetTagCollection> bTagCollectionSSVHE;
    iEvent.getByLabel("simpleSecondaryVertexHighEffBJetTags", bTagCollectionSSVHE);
    const reco::JetTagCollection & bTagsSSVHE = *(bTagCollectionSSVHE.product());

    edm::Handle<reco::JetTagCollection> bTagCollectionJBP;
    iEvent.getByLabel("jetBProbabilityBJetTags", bTagCollectionJBP);
    const reco::JetTagCollection & bTagsJBP = *(bTagCollectionJBP.product());

    edm::Handle<reco::JetTagCollection> bTagCollectionCSV;
    iEvent.getByLabel("combinedSecondaryVertexBJetTags", bTagCollectionCSV);
    const reco::JetTagCollection & bTagsCSV = *(bTagCollectionCSV.product());

    edm::Handle<reco::JetTagCollection> bTagCollectionCSVMVA;
    iEvent.getByLabel("combinedSecondaryVertexBJetTags", bTagCollectionCSVMVA);
    const reco::JetTagCollection & bTagsCSVMVA = *(bTagCollectionCSVMVA.product());

    // Jet flavor matching //
    Handle<JetFlavourMatchingCollection> bJetFlavourMC;
    if (!isRealData) {
        iEvent.getByLabel("JetbyValAlgo", bJetFlavourMC);
    }

    Handle<vector<reco::PFJet> > jets;
    iEvent.getByLabel(jetTag_, jets);

    for (vector<reco::PFJet>::const_iterator iJet = jets->begin(); iJet != jets->end(); ++iJet) {

        if (iJet->pt() < 10.) continue;

        TCJet* jetCon (new ((*recoJets)[jetCount]) TCJet);

        jetCon->SetPxPyPzE(iJet->px(), iJet->py(), iJet->pz(), iJet->energy());
        jetCon->SetVtx(0., 0., 0.);
        jetCon->SetChHadFrac(iJet->chargedHadronEnergyFraction());
        jetCon->SetNeuHadFrac(iJet->neutralHadronEnergyFraction());
        jetCon->SetChEmFrac(iJet->chargedEmEnergyFraction());
        jetCon->SetNeuEmFrac(iJet->neutralEmEnergyFraction());
        jetCon->SetNumConstit(iJet->chargedMultiplicity() + iJet->neutralMultiplicity());
        jetCon->SetNumChPart(iJet->chargedMultiplicity());

        //jetCon->SetJetFlavor(iJet->partonFlavour());
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
        jetCon->SetBDiscriminatorMap("CSVMVA", MatchBTagsToJets(bTagsCSVMVA, *iJet));

        /////////////////////
        // Get Hgg Id vars //
        /////////////////////

        PileupJetIdentifier puIdentifier;
        // giving uncorrected input, must double check on this
        float jec = 1.;
        // jet corrector
        if( jecCor.get() == 0 ) {
            initJetEnergyCorrector( iSetup, iEvent.isRealData() );
        }

        jecCor->setJetPt(iJet->pt());
        jecCor->setJetEta(iJet->eta());
        jecCor->setJetA(iJet->jetArea());
        jecCor->setRho(rhoFactor);
        jec = jecCor->getCorrection();
        //cout<<"jec:\t"<<jec<<endl;
        VertexCollection::const_iterator vtx;
        const VertexCollection & vertexes = *(primaryVtcs.product());
        vtx = vertexes.begin();
        while( vtx != vertexes.end() && ( vtx->isFake() || vtx->ndof() < 4 ) ) {
            ++vtx;
        }
        if( vtx == vertexes.end() ) { vtx = vertexes.begin(); }
        puIdentifier = myPUJetID->computeIdVariables(&(*iJet), jec,  &(*vtx),  vertexes, false);
        //cout<<"betaStarClassic:\t"<<puIdentifier.betaStarClassic()<<"\t"<<"dR2Mean:\t"<<puIdentifier.dR2Mean()<<endl;
        jetCon->SetBetaStarClassic(puIdentifier.betaStarClassic());
        jetCon->SetDR2Mean(puIdentifier.dR2Mean());

        ++jetCount;
    }


    /////////////
    // Get MET //
    /////////////


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

    edm::Handle< edm::View<reco::PFMET> > pfMEThandle;
    iEvent.getByLabel("pfMet", pfMEThandle);
    //Significance
    float significance = (pfMEThandle->front()).significance();
    float sigmaX2 = (pfMEThandle->front()).getSignificanceMatrix()(0,0);
    recoMET->SetSignificance( significance );
    recoMET->SetSigmaX2( sigmaX2 );


    ///////////////
    // Get muons //
    ///////////////



    Handle<vector<reco::Muon> > muons;
    iEvent.getByLabel(muonTag_, muons);

    for (vector<reco::Muon>::const_iterator iMuon = muons->begin(); iMuon != muons->end(); ++iMuon) {
        //if (!iMuon->isGlobalMuon() || iMuon->pt() < 3.) continue;
        //if (iMuon->pt() < 3.) continue;
        //if (!iMuon->isGlobalMuon()) continue;

        TCMuon* myMuon = new ((*recoMuons)[muCount]) TCMuon;

        myMuon->SetPxPyPzE(iMuon->px(), iMuon->py(), iMuon->pz(), iMuon->energy());
        myMuon->SetCharge(iMuon->charge());

        myMuon->SetPF(iMuon->isPFMuon());
        myMuon->SetIsGLB(iMuon->isGlobalMuon());
        myMuon->SetIsTRK(iMuon->isTrackerMuon());


        myMuon->SetIsGood(muon::isGoodMuon(*iMuon, muon::TMOneStationTight));
        myMuon->SetIsGoodLoose(muon::isGoodMuon(*iMuon, muon::TMOneStationLoose));

        myMuon->SetIsArbitrated(muon::isGoodMuon(*iMuon, muon::AllArbitrated));
        myMuon->SetIsTrkArbitrated(muon::isGoodMuon(*iMuon, muon::TrackerMuonArbitrated));


        if (primaryVtcs->size()>0){
            myMuon->SetIsTight(muon::isTightMuon(*iMuon, *primaryVtcs->begin()));
            myMuon->SetIsSoft( muon::isSoftMuon( *iMuon, *primaryVtcs->begin()));
        }
        else{
            myMuon->SetIsTight(0);
            myMuon->SetIsSoft(0);
        }

        myMuon->SetCaloComp(iMuon->caloCompatibility());
        myMuon->SetSegComp(muon::segmentCompatibility(*iMuon));
        myMuon->SetNumberOfMatchedStations(iMuon->numberOfMatchedStations());
        myMuon->SetNumberOfMatches(iMuon->numberOfMatches());

        if (iMuon->isGlobalMuon()){
            myMuon->SetNormalizedChi2(       iMuon->globalTrack()->normalizedChi2());
            myMuon->SetNumberOfValidMuonHits(iMuon->globalTrack()->hitPattern().numberOfValidMuonHits());
        }
        else{
            myMuon->SetNormalizedChi2(-1);
            myMuon->SetNumberOfValidMuonHits(-1);
        }

        if (iMuon->isTrackerMuon()){
            myMuon->SetVtx(iMuon->track()->vx(),iMuon->track()->vy(),iMuon->track()->vz());
            myMuon->SetPtError(iMuon->track()->ptError());

            myMuon->SetTrackLayersWithMeasurement(iMuon->track()->hitPattern().trackerLayersWithMeasurement());
            myMuon->SetPixelLayersWithMeasurement(iMuon->innerTrack()->hitPattern().pixelLayersWithMeasurement());
            myMuon->SetNumberOfValidPixelHits(    iMuon->innerTrack()->hitPattern().numberOfValidPixelHits());
            myMuon->SetNormalizedChi2_tracker(    iMuon->innerTrack()->normalizedChi2());
            myMuon->SetNumberOfValidTrackerHits(iMuon->track()->hitPattern().numberOfValidTrackerHits());
            myMuon->SetNumberOfLostPixelHits(   iMuon->track()->hitPattern().numberOfLostPixelHits());
            myMuon->SetNumberOfLostTrackerHits( iMuon->track()->hitPattern().numberOfLostTrackerHits());
        }
        else{
            myMuon->SetVtx(-1,-1,-1);
            myMuon->SetPtError(-1);
            myMuon->SetTrackLayersWithMeasurement(-1);
            myMuon->SetNumberOfValidPixelHits(-1);
            myMuon->SetNormalizedChi2_tracker(-1);
            myMuon->SetNumberOfValidTrackerHits(-1);
            myMuon->SetNumberOfLostPixelHits(-1);
            myMuon->SetNumberOfLostTrackerHits(-1);
        }
        // Set isolation map values
        // Detector-based isolation

        /*
           myMuon->SetIsoMap("NTracks_R03", iMuon->isolationR03().nTracks);
           myMuon->SetIsoMap("EmIso_R03",   iMuon->isolationR03().emEt);
           myMuon->SetIsoMap("HadIso_R03",  iMuon->isolationR03().hadEt);
           myMuon->SetIsoMap("SumPt_R03",   iMuon->isolationR03().sumPt);

           myMuon->SetIsoMap("NTracks_R05", iMuon->isolationR05().nTracks);
           myMuon->SetIsoMap("EmIso_R05",   iMuon->isolationR05().emEt);
           myMuon->SetIsoMap("HadIso_R05",  iMuon->isolationR05().hadEt);
           myMuon->SetIsoMap("SumPt_R05",   iMuon->isolationR05().sumPt);
         */
        // PF-based isolation
        myMuon->SetIdMap("Iso_pfPUPt_R03",      iMuon->pfIsolationR03().sumPUPt);
        myMuon->SetIdMap("Iso_pfPhotonEt_R03",  iMuon->pfIsolationR03().sumPhotonEt);
        myMuon->SetIdMap("Iso_pfChargedPt_R03", iMuon->pfIsolationR03().sumChargedParticlePt);
        myMuon->SetIdMap("Iso_pfChargedHadronPt_R03", iMuon->pfIsolationR03().sumChargedHadronPt);
        myMuon->SetIdMap("Iso_pfNeutralHadronEt_R03", iMuon->pfIsolationR03().sumNeutralHadronEt);
        /*
           myMuon->SetIsoMap("pfPUPt_R04",      iMuon->pfIsolationR04().sumPUPt);
           myMuon->SetIsoMap("pfPhotonEt_R04",  iMuon->pfIsolationR04().sumPhotonEt);
           myMuon->SetIsoMap("pfChargedPt_R04", iMuon->pfIsolationR04().sumChargedParticlePt);
           myMuon->SetIsoMap("pfChargedHadronPt_R04", iMuon->pfIsolationR04().sumChargedHadronPt);
           myMuon->SetIsoMap("pfNeutralHadronEt_R04", iMuon->pfIsolationR04().sumNeutralHadronEt);
         */

        myMuon->SetPfIsoPU(iMuon->pfIsolationR04().sumPUPt);
        myMuon->SetPfIsoChargedPart(iMuon->pfIsolationR04().sumChargedParticlePt);
        myMuon->SetPfIsoChargedHad( iMuon->pfIsolationR04().sumChargedHadronPt);
        myMuon->SetPfIsoNeutral(iMuon->pfIsolationR04().sumNeutralHadronEt);
        myMuon->SetPfIsoPhoton( iMuon->pfIsolationR04().sumPhotonEt);

        //cout << "\t" << iMuon->pt() << ", " << iMuon->eta() << endl;
        // Match muon to trigger object //
        if (saveTriggerObj_){
          for (unsigned j = 0; j < triggerObjects.size(); ++j) {
              float deltaR = triggerObjects[j].DeltaR(*myMuon);
              if (deltaR < 0.3 && fabs(triggerObjects[j].GetId()) == 13) {
                  myMuon->SetTriggered(true);
                  break;
              } 
          }
        }
        muCount++;
        //cout<<*myMuon<<endl;
    }


    ///////////////////
    // Get electrons //
    ///////////////////


    edm::Handle<reco::ConversionCollection> hConversions;
    iEvent.getByLabel("allConversions", hConversions);


    Handle<reco::GsfElectronCollection > electrons;
    iEvent.getByLabel(electronTag_, electrons);

    Handle<reco::GsfElectronCollection > calibratedElectrons;
    iEvent.getByLabel(edm::InputTag("calibratedElectrons","calibratedGsfElectrons"), calibratedElectrons);

    edm::Handle<edm::ValueMap<float>> mvaTrigV0_handle;
    iEvent.getByLabel("mvaTrigV0", mvaTrigV0_handle);
    const edm::ValueMap<float> ele_mvaTrigV0 = (*mvaTrigV0_handle.product());

    edm::Handle<edm::ValueMap<float>> mvaNonTrigV0_handle;
    iEvent.getByLabel("mvaNonTrigV0", mvaNonTrigV0_handle);
    const edm::ValueMap<float> ele_mvaNonTrigV0 = (*mvaNonTrigV0_handle.product());

    edm::Handle<edm::ValueMap<double>> regEne_handle;
    iEvent.getByLabel(edm::InputTag("eleRegressionEnergy","eneRegForGsfEle"), regEne_handle);
    const edm::ValueMap<double> ele_regEne = (*regEne_handle.product());

    edm::Handle<edm::ValueMap<double>> regErr_handle;
    iEvent.getByLabel(edm::InputTag("eleRegressionEnergy","eneErrorRegForGsfEle"), regErr_handle);
    const edm::ValueMap<double> ele_regErr = (*regErr_handle.product());



    // DO NOT delete it yet, it maybe useful later for Dalitz electron isolation
    //This stuff is for modified isolation for close electrons,
    //following prescription here: https://twiki.cern.ch/twiki/bin/viewauth/CMS/BoostedZToEEModIso
    edm::Handle<edm::ValueMap<double> > h_modElectronIso_Tk;
    edm::Handle<edm::ValueMap<double> > h_modElectronIso_Ecal;
    edm::Handle<edm::ValueMap<double> > h_modElectronIso_HcalD1;
    iEvent.getByLabel("modElectronIso","track",      h_modElectronIso_Tk);
    iEvent.getByLabel("modElectronIso","ecal",       h_modElectronIso_Ecal);
    iEvent.getByLabel("modElectronIso","hcalDepth1", h_modElectronIso_HcalD1);
    const edm::ValueMap<double> modElectronIso_Tk     = (*h_modElectronIso_Tk.product());
    const edm::ValueMap<double> modElectronIso_Ecal   = (*h_modElectronIso_Ecal.product());
    const edm::ValueMap<double> modElectronIso_HcalD1 = (*h_modElectronIso_HcalD1.product());


    // These are the Value Maps for Iso deposits:

    edm::Handle<edm::ValueMap<double>> iso_charged_handle;
    iEvent.getByLabel(edm::InputTag("elPFIsoValueCharged04PFIdPFIso",""), iso_charged_handle);
    const edm::ValueMap<double> ele_iso_charged = (*iso_charged_handle.product());

    edm::Handle<edm::ValueMap<double>> iso_gamma_handle;
    iEvent.getByLabel(edm::InputTag("elPFIsoValueGamma04PFIdPFIso",""), iso_gamma_handle);
    const edm::ValueMap<double> ele_iso_gamma = (*iso_gamma_handle.product());

    edm::Handle<edm::ValueMap<double>> iso_neutral_handle;
    iEvent.getByLabel(edm::InputTag("elPFIsoValueNeutral04PFIdPFIso",""), iso_neutral_handle);
    const edm::ValueMap<double> ele_iso_neutral = (*iso_neutral_handle.product());

    Int_t eee=0;
    for (vector<reco::GsfElectron>::const_iterator iElectron = electrons->begin(); iElectron != electrons->end(); ++iElectron) {
        eee++;
        if (iElectron->pt() < 5) continue;

        TCElectron* myElectron = new ((*recoElectrons)[eleCount]) TCElectron;

        // Basic physics object info
        myElectron->SetPtEtaPhiM(iElectron->pt(), iElectron->eta(), iElectron->phi(), iElectron->mass());
        myElectron->SetVtx(iElectron->gsfTrack()->vx(),iElectron->gsfTrack()->vy(),iElectron->gsfTrack()->vz());
        myElectron->SetCharge(iElectron->charge());

        // Fiducial variables
        myElectron->SetIsEB(iElectron->isEB());
        myElectron->SetIsEE(iElectron->isEE());
        myElectron->SetIsInGap(iElectron->isGap());


        // Electron ID variables

        //Methods that are availabel for the electrons can be found here:
        //http://cmslxr.fnal.gov/lxr/source/DataFormats/EgammaCandidates/interface/GsfElectron.h?v=CMSSW_5_3_11

        myElectron->SetR9(     iElectron->r9());
        myElectron->SetFBrem(  iElectron->fbrem());
        myElectron->SetEoP(    iElectron->eSuperClusterOverP());
        myElectron->SetEoPout( iElectron->eEleClusterOverPout());
        myElectron->SetESeedOverP(iElectron->eSeedClusterOverP());

        // ***** >> Switching to officially recommended method:  <<<<<<
        myElectron->SetHadOverEm(iElectron->hcalOverEcalBc());
        // !!!!!!!!!
        // QUESTION: Does the eleIsolator below returns the recommended isolation for Hcal??
        //!!!!!!!!!!

        //cout<<"H/E compare: hcalOverEcalBc = "<<iElectron->hcalOverEcalBc()<<"   hadronicOverEm = "<<iElectron->hadronicOverEm()<<endl;

        myElectron->SetSCEta(  iElectron->superCluster()->eta());
        myElectron->SetSCPhi(  iElectron->superCluster()->phi());

        //*** --> Notice that previously some variables were defined in the IdMap, differently:
        //one has to perform a selections on analysis level to recover this behaviour:
        //(cut at over/underflow values)

        myElectron->SetSCDeltaEta(   iElectron->deltaEtaSuperClusterTrackAtVtx());
        myElectron->SetSCDeltaPhi(   iElectron->deltaPhiSuperClusterTrackAtVtx());
        myElectron->SetSigmaIEtaIEta(iElectron->sigmaIetaIeta());
        myElectron->SetSigmaIPhiIPhi(iElectron->sigmaIphiIphi());
        myElectron->SetSCEtaWidth(   iElectron->superCluster()->etaWidth());
        myElectron->SetSCPhiWidth(   iElectron->superCluster()->phiWidth());
        myElectron->SetSCEnergy(     iElectron->superCluster()->energy());
        if (iElectron->superCluster()->rawEnergy()!=0)
            myElectron->SetPreShowerOverRaw(iElectron->superCluster()->preshowerEnergy() / iElectron->superCluster()->rawEnergy());


        myElectron->SetE1x5(iElectron->e1x5());
        myElectron->SetE2x5(iElectron->e2x5Max());
        myElectron->SetE5x5(iElectron->e5x5());

        myElectron->SetDeltaEtaSeedCluster(iElectron->deltaEtaSeedClusterTrackAtCalo());
        myElectron->SetDeltaPhiSeedCluster(iElectron->deltaPhiSeedClusterTrackAtCalo());

        //std::vector vCov = iElectron->superCluster()->localCovariances(*(iElectron->superCluster()->seed())) ;


        myElectron->SetEoP(iElectron->eSuperClusterOverP());

        // ** *************
        // Assosited GSF tracks:
        // ** ************

        TCElectron::Track *t = new TCElectron::Track();
        t->SetXYZM(iElectron->gsfTrack()->px(), iElectron->gsfTrack()->py(), iElectron->gsfTrack()->pz(),  0);
        t->SetVtx(iElectron->gsfTrack()->vx(),  iElectron->gsfTrack()->vy(), iElectron->gsfTrack()->vz());
        t->SetCharge(iElectron->gsfTrack()->chargeMode());
        t->SetNormalizedChi2(iElectron->gsfTrack()->normalizedChi2());
        t->SetPtError(iElectron->gsfTrack()->ptError());
        TCTrack::ConversionInfo convInfo = ntupleProducer::CheckForConversions(hConversions, iElectron->gsfTrack(),
                vertexBeamSpot.position(), (*myVtxRef).position());
        t->SetConversionInfo(convInfo);
        //This is the main track, directly assosiated with an Electron
        myElectron->AddTrack(*t);

        myElectron->SetPtError(iElectron->gsfTrack()->ptError());

        //Int_t ntr=0;
        //Adding more tracks from the ambiguos collection:
        for (reco::GsfTrackRefVector::const_iterator gtr = iElectron->ambiguousGsfTracksBegin(); gtr != iElectron->ambiguousGsfTracksEnd(); ++gtr)
        {
            //ntr++;
            //cout<<ntr<<" ambigious loop pt="<<(*gtr)->pt()<<" eta="<<(*gtr)->eta()<<" "<<" phi="<<(*gtr)->phi()<<endl;

            t->SetXYZM((*gtr)->px(), (*gtr)->py(), (*gtr)->pz(),  0);
            t->SetVtx((*gtr)->vx(), (*gtr)->vy(), (*gtr)->vz());
            t->SetCharge((*gtr)->chargeMode());
            t->SetNormalizedChi2((*gtr)->normalizedChi2());
            t->SetPtError((*gtr)->ptError());
            //re-using the same object
            convInfo = ntupleProducer::CheckForConversions(hConversions, *gtr,
                    vertexBeamSpot.position(), (*myVtxRef).position());
            t->SetConversionInfo(convInfo);
            myElectron->AddTrack(*t);

        }

        myElectron->SetConversionDcot(iElectron->convDcot());
        myElectron->SetConversionDist(iElectron->convDist());
        myElectron->SetConversionRadius(iElectron->convRadius());

        bool validKF= false;
        reco::TrackRef myTrackRef = iElectron->closestCtfTrackRef();
        validKF = (myTrackRef.isAvailable());
        validKF = (myTrackRef.isNonnull());

        if (validKF){
            myElectron->SetTrackerLayersWithMeasurement( myTrackRef->hitPattern().trackerLayersWithMeasurement());
            myElectron->SetNormalizedChi2Kf( myTrackRef->normalizedChi2());
            myElectron->SetNumberOfValidHits(myTrackRef->numberOfValidHits());

        }
        else{
            myElectron->SetTrackerLayersWithMeasurement(-1);
            myElectron->SetNormalizedChi2Kf(-1);
            myElectron->SetNumberOfValidHits(-1);
        }


        Handle<reco::VertexCollection> thePrimaryVertexColl;
        iEvent.getByLabel("offlinePrimaryVertices",thePrimaryVertexColl);

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

        float ip3d    = -999.0;
        float ip3derr = 1.0;
        float ip3dSig = 0.0;

        edm::ESHandle<TransientTrackBuilder> builder;
        iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", builder);
        TransientTrackBuilder thebuilder = *(builder.product());

        if (iElectron->gsfTrack().isNonnull()) {
            const double gsfsign   = ( (-iElectron->gsfTrack()->dxy(pv->position()))   >=0 ) ? 1. : -1.;

            const reco::TransientTrack &tt = thebuilder.build(iElectron->gsfTrack());
            const std::pair<bool,Measurement1D> &ip3dpv =  IPTools::absoluteImpactParameter3D(tt,*pv);
            if (ip3dpv.first) {
                ip3d = gsfsign*ip3dpv.second.value();
                ip3derr = ip3dpv.second.error();
                ip3dSig = ip3d/ip3derr;
            }
        }

        myElectron->SetIP3d(ip3d);
        myElectron->SetIP3dSig(ip3dSig);

        myElectron->SetNumberOfValidPixelHits(  iElectron->gsfTrack()->hitPattern().numberOfValidPixelHits());
        myElectron->SetNumberOfValidTrackerHits(iElectron->gsfTrack()->hitPattern().numberOfValidTrackerHits());
        myElectron->SetNumberOfLostPixelHits(   iElectron->gsfTrack()->hitPattern().numberOfLostPixelHits());
        myElectron->SetNumberOfLostTrackerHits( iElectron->gsfTrack()->hitPattern().numberOfLostTrackerHits());

        myElectron->SetInverseEnergyMomentumDiff(fabs((1/iElectron->ecalEnergy()) - (1/iElectron->trackMomentumAtVtx().R())));

        // Conversion information
        // See definition from here: https://twiki.cern.ch/twiki/bin/view/CMS/ConversionTools
        bool passConvVeto = !(ConversionTools::hasMatchedConversion(*iElectron,hConversions,vertexBeamSpot.position()));
        myElectron->SetPassConversionVeto(passConvVeto);
        myElectron->SetConversionMissHits(iElectron->gsfTrack()->trackerExpectedHitsInner().numberOfHits());

        //Variables for Triggering MVA pre-selection
        myElectron->SetIdMap("hadronicOverEm",      iElectron->hadronicOverEm());
        myElectron->SetIdMap("dr03TkSumPt",         iElectron->dr03TkSumPt());
        myElectron->SetIdMap("dr03EcalRecHitSumEt", iElectron->dr03EcalRecHitSumEt());
        myElectron->SetIdMap("dr03HcalTowerSumEt",  iElectron->dr03HcalTowerSumEt());
        myElectron->SetIdMap("gsf_numberOfLostHits",iElectron->gsfTrack()->trackerExpectedHitsInner().numberOfLostHits());

        // Effective area for rho PU corrections 
        float AEff03 = ElectronEffectiveArea::GetElectronEffectiveArea(ElectronEffectiveArea::kEleGammaAndNeutralHadronIso03, iElectron->eta(), ElectronEffectiveArea::kEleEAData2012);
        float AEff04 = ElectronEffectiveArea::GetElectronEffectiveArea(ElectronEffectiveArea::kEleGammaAndNeutralHadronIso04, iElectron->eta(), ElectronEffectiveArea::kEleEAData2012);
        myElectron->SetIdMap("EffArea_R03", AEff03);
        myElectron->SetIdMap("EffArea_R04", AEff04);
        myElectron->SetEffArea(AEff04);


        //MVA output:
        float m_old = ele_mvaTrigV0.get(eee-1);
        myElectron->SetMvaID_Old(m_old);
        float m_HZZ = ele_mvaNonTrigV0.get(eee-1);
        myElectron->SetMvaID_HZZ(m_HZZ);

        //cout<<eee<<"  mva0 = "<<m<<endl;

        //Regression energy
        double ene = ele_regEne.get(eee-1);
        double err = ele_regErr.get(eee-1);
        myElectron->SetEnergyRegression(ene);
        myElectron->SetEnergyRegressionErr(err);

        myElectron->SetIdMap("modElectronIso_Tk",    modElectronIso_Tk.get(eee-1));
        myElectron->SetIdMap("modElectronIso_Ecal",  modElectronIso_Ecal.get(eee-1));
        myElectron->SetIdMap("modElectronIso_HcalD1",modElectronIso_HcalD1.get(eee-1));

        //cout<<eee-1<<"Bosted iso vars: tk="<<modElectronIso_Tk.get(eee-1)
        //  <<"   ecal="<<modElectronIso_Ecal.get(eee-1)
        //  <<" hcald1="<<modElectronIso_HcalD1.get(eee-1)<<endl;

        const reco::GsfElectron &iElectronTmp((*calibratedElectrons)[eee-1]);

        //cout<<"ielectron , pt ="<<iElectron->pt()<<" eta="<<iElectron->eta()<<endl;
        //cout<<"  calibra , pt ="<<iElectronTmp.pt()<<" eta="<<iElectronTmp.eta()<<endl;

        TLorentzVector tmpP4;
        tmpP4.SetPtEtaPhiE(iElectronTmp.pt(), iElectronTmp.eta(), iElectronTmp.phi(), iElectronTmp.energy());
        myElectron->SetRegressionMomCombP4(tmpP4);

        //3 types of Pf Iso

        // 1 PFIsolationEstimator
        eleIsolator.fGetIsolation(&(*iElectron), &thePfColl, myVtxRef, primaryVtcs);

        // 2 CommonTools/ParticleFlow
        double iso_charged = ele_iso_charged.get(eee-1);
        double iso_neutral = ele_iso_neutral.get(eee-1);
        double iso_gamma   = ele_iso_gamma.get(eee-1);

        // 3 SCFootprint removal
        edm::ParameterSet myConfig;
        myConfig.addUntrackedParameter("isolation_cone_size_forSCremoval",SCFPRemovalConeEl_);
        SuperClusterFootprintRemoval remover(iEvent,iSetup,myConfig);
        PFIsolation_struct mySCFPstruct = remover.PFIsolation(iElectron->superCluster(),edm::Ptr<Vertex>(primaryVtcs, 0));

        //cout<<"Electron isolation printout for electron # "<<eee-1<<"   pt ="<<iElectron->pt()<<endl;
        //cout<<"old charged = "<<eleIsolator.getIsolationCharged()<<"\t gamma = "<<eleIsolator.getIsolationPhoton()<<"\t neutral = "<<eleIsolator.getIsolationNeutral()<<endl;
        //cout<<"new charged = "<<iso_charged<<"\t gamma = "<<iso_gamma<<"\t neutral = "<<iso_neutral<<endl;
        //cout<<"alt charged = "<<mySCFPstruct.chargediso<<"\t gamma = "<<mySCFPstruct.photoniso<<"\t neutral = "<<mySCFPstruct.neutraliso<<endl;

        myElectron->SetPfIsoCharged(eleIsolator.getIsolationCharged());
        myElectron->SetPfIsoNeutral(eleIsolator.getIsolationNeutral());
        myElectron->SetPfIsoPhoton(eleIsolator.getIsolationPhoton());

        myElectron->SetIdMap("chIso_std",  iso_charged);
        myElectron->SetIdMap("neuIso_std", iso_neutral);
        myElectron->SetIdMap("phoIso_std", iso_gamma);

        myElectron->SetIdMap("chIso_scfp", mySCFPstruct.chargediso);
        myElectron->SetIdMap("neuIso_scfp", mySCFPstruct.neutraliso);
        myElectron->SetIdMap("phoIso_scfp", mySCFPstruct.photoniso);

        // Match muon to trigger object //
        if(saveTriggerObj_){
          for (unsigned j = 0; j < triggerObjects.size(); ++j) {
              float deltaR    = triggerObjects[j].DeltaR(*myElectron);
              //float deltaPt   = fabs(triggerObjects[j].Pt() - myElectron->Pt())/myElectron->Pt();
              if (deltaR < 0.1 && fabs(triggerObjects[j].GetId()) == 11) {
                  //cout << triggerObjects[j].GetHLTName() << "\t" << triggerObjects[j].GetModuleName() << "\t" << deltaR << "\t" << deltaPt << endl;;
                  myElectron->SetTriggered(true);
                  break;
              } 
          }
        }
        eleCount++;
    }


    /////////////////
    // Get photons //
    /////////////////

    // ES geometry
    ESHandle<CaloGeometry> geoHandle;
    iSetup.get<CaloGeometryRecord>().get(geoHandle);
    const CaloSubdetectorGeometry *geometry = geoHandle->getSubdetectorGeometry(DetId::Ecal, EcalPreshower);
    const CaloSubdetectorGeometry *& geometry_p = geometry;

    if (geometry) topology_p.reset(new EcalPreshowerTopology(geoHandle));

    // make the map of rechits
    Handle<EcalRecHitCollection> ESRecHits;
    iEvent.getByLabel(esReducedRecHitCollection_,ESRecHits);

    rechits_map_.clear();
    if (ESRecHits.isValid()) {
        EcalRecHitCollection::const_iterator it;
        for (it = ESRecHits->begin(); it != ESRecHits->end(); ++it) {
            // remove bad ES rechits
            if (it->recoFlag()==1 || it->recoFlag()==14 || (it->recoFlag()<=10 && it->recoFlag()>=5)) continue;
            //Make the map of DetID, EcalRecHit pairs
            rechits_map_.insert(std::make_pair(it->id(), *it));
        }
    }

    Handle<EcalRecHitCollection> Brechit;
    iEvent.getByLabel("reducedEcalRecHitsEB",Brechit);
    //const EcalRecHitCollection* barrelRecHits= Brechit.product();

    Handle<vector<reco::Photon> > photons;
    iEvent.getByLabel(photonTag_, photons);

    edm::Handle<reco::GsfElectronCollection> hElectrons;
    iEvent.getByLabel("gsfElectrons", hElectrons);


    for (vector<reco::Photon>::const_iterator iPhoton = photons->begin(); iPhoton != photons->end() ; ++iPhoton) {
        TCPhoton* myPhoton = new ((*recoPhotons)[photonCount]) TCPhoton();

        size_t rightRecoPho = -1;
        for (size_t iv = 0; iv < photons->size(); ++iv) {
            reco::PhotonRef recophoRef2(photons, iv);
            if (deltaR(iPhoton->eta(), iPhoton->phi(), recophoRef2->eta(), recophoRef2->phi()) < 0.01) rightRecoPho = iv;
        }
        reco::PhotonRef recophoRef(photons, rightRecoPho);

        if (savePhoCrystals_)
        {
            //Crystal Info:
            std::vector< std::pair<DetId, float> >  PhotonHit_DetIds  = iPhoton->superCluster()->hitsAndFractions();
            std::vector<TCEGamma::CrystalInfo> crystalinfo_container;
            crystalinfo_container.clear();
            TCPhoton::CrystalInfo crystal = {};
            float timing_avg =0.0;
            int ncrys   = 0;
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

            myPhoton->SetNCrystals(crystalinfo_container.size());

            for (unsigned int y =0; y < crystalinfo_container.size() && y < 100;y++){
                myPhoton->AddCrystal(crystalinfo_container[y]);
            }

            /*
               vector<TCPhoton::CrystalInfo> savedCrystals = myPhoton->GetCrystalVect();
               for (int y = 0; y< myPhoton->GetNCrystals();y++){
               std::cout << "savedCrystals[y].time : " << savedCrystals[y].time << std::endl;
               std::cout << "savedCrystals[y].timeErr : " << savedCrystals[y].timeErr << std::endl;
               std::cout << "savedCrystals[y].energy : " << savedCrystals[y].energy <<std::endl;
               std::cout << "savedCrystals[y].ieta: " << savedCrystals[y].ieta << std::endl;

               std::cout << "savedCrystals[y].rawId: " << savedCrystals[y].rawId <<std::endl;
               }
             */

        }

        myPhoton->SetPxPyPzE(iPhoton->px(), iPhoton->py(), iPhoton->pz(), iPhoton->p());
        myPhoton->SetVtx(iPhoton->vx(), iPhoton->vy(), iPhoton->vz());

        // ID variables
        //Methods that are availabel for the electrons can be found here:
        //http://cmslxr.fnal.gov/lxr/source/DataFormats/EgammaCandidates/interface/Photon.h?v=CMSSW_5_3_11

        vector<float> phoCov;
        const reco::CaloClusterPtr phoSeed = iPhoton->superCluster()->seed();
        phoCov = lazyTool->localCovariances(*phoSeed);

        myPhoton->SetHadOverEm(iPhoton->hadTowOverEm());
        myPhoton->SetR9(iPhoton->r9());
        myPhoton->SetTrackVeto(iPhoton->hasPixelSeed());

        myPhoton->SetSCEta(iPhoton->superCluster()->eta());
        myPhoton->SetSCPhi(iPhoton->superCluster()->phi());
        myPhoton->SetSigmaIEtaIEta(iPhoton->sigmaIetaIeta());
        myPhoton->SetSigmaIEtaIPhi(phoCov[1]);
        myPhoton->SetSigmaIPhiIPhi(phoCov[2]);

        myPhoton->SetSCEtaWidth(  iPhoton->superCluster()->etaWidth());
        myPhoton->SetSCPhiWidth(  iPhoton->superCluster()->phiWidth());

        myPhoton->SetSCEnergy(iPhoton->superCluster()->energy());
        myPhoton->SetSCRawEnergy(iPhoton->superCluster()->rawEnergy());
        myPhoton->SetSCPSEnergy(iPhoton->superCluster()->preshowerEnergy());

        if (iPhoton->superCluster()->rawEnergy()!=0)
            myPhoton->SetPreShowerOverRaw(iPhoton->superCluster()->preshowerEnergy() / iPhoton->superCluster()->rawEnergy());


        myPhoton->SetE1x3(lazyTool->e1x3(*phoSeed));
        myPhoton->SetE1x5(iPhoton->e1x5());
        myPhoton->SetE2x2(lazyTool->e2x2(*phoSeed));
        myPhoton->SetE2x5(iPhoton->e2x5()); // <<-
        // How come these two aren't  the same?!
        myPhoton->SetE2x5Max(lazyTool->e2x5Max(*phoSeed)); //<<-
        //if (iPhoton->e2x5() != lazyTool->e2x5Max(*phoSeed))
        //cout<<"No, it's not the same! : "<< iPhoton->e2x5()<<" vs "<< lazyTool->e2x5Max(*phoSeed)<<endl;

        myPhoton->SetE5x5(iPhoton->e5x5());

        // PF Iso for photons
        phoIsolator.fGetIsolation(&(*iPhoton),&thePfColl, myVtxRef, primaryVtcs);
        myPhoton->SetPfIsoCharged(phoIsolator.getIsolationCharged());
        myPhoton->SetPfIsoNeutral(phoIsolator.getIsolationNeutral());
        myPhoton->SetPfIsoPhoton( phoIsolator.getIsolationPhoton());

        // CiC track Iso
        vector<float> vtxIsolations02 = cicPhotonId_->pfTkIsoWithVertex(recophoRef, 0.2, 0.02, 0.02, 0.0, 0.2, 0.1, reco::PFCandidate::h);
        vector<float> vtxIsolations03 = cicPhotonId_->pfTkIsoWithVertex(recophoRef, 0.3, 0.02, 0.02, 0.0, 0.2, 0.1, reco::PFCandidate::h);
        vector<float> vtxIsolations04 = cicPhotonId_->pfTkIsoWithVertex(recophoRef, 0.4, 0.02, 0.02, 0.0, 0.2, 0.1, reco::PFCandidate::h);

        myPhoton->SetCiCPF4chgpfIso02(vtxIsolations02);
        myPhoton->SetCiCPF4chgpfIso03(vtxIsolations03);
        myPhoton->SetCiCPF4chgpfIso04(vtxIsolations04);

        /*
           vector<float> testIso = myPhoton->CiCPF4chgpfIso02();
           cout<<"test CiC: "<<testIso[0]<<endl;
         */


        // detector-based isolation
        myPhoton->SetIdMap("EmIso_R03",  (iPhoton->ecalRecHitSumEtConeDR03()));
        myPhoton->SetIdMap("HadIso_R03", (iPhoton->hcalTowerSumEtConeDR03()));
        myPhoton->SetIdMap("TrkIso_R03", (iPhoton->trkSumPtHollowConeDR03()));


        myPhoton->SetIdMap("EmIso_R04",  (iPhoton->ecalRecHitSumEtConeDR04()));
        myPhoton->SetIdMap("HadIso_R04", (iPhoton->hcalTowerSumEtConeDR04()));
        myPhoton->SetIdMap("TrkIso_R04", (iPhoton->trkSumPtHollowConeDR04()));


        // Hcal isolation for 2012
        //myPhoton->SetIdMap("HadIso_R03",iPhoton->hcalTowerSumEtConeDR03() +
        //        (iPhoton->hadronicOverEm() - iPhoton->hadTowOverEm())*iPhoton->superCluster()->energy()/cosh(iPhoton->superCluster()->eta()));
        //myPhoton->SetIdMap("HadIso_R04",iPhoton->hcalTowerSumEtConeDR04() +
        //        (iPhoton->hadronicOverEm() - iPhoton->hadTowOverEm())*iPhoton->superCluster()->energy()/cosh(iPhoton->superCluster()->eta()));

        //Footprint removal
        edm::ParameterSet myConfig;
        myConfig.addUntrackedParameter("isolation_cone_size_forSCremoval",SCFPRemovalConePho_);
        SuperClusterFootprintRemoval remover(iEvent,iSetup,myConfig);
        PFIsolation_struct mySCFPstruct = remover.PFIsolation(iPhoton->superCluster(),edm::Ptr<Vertex>(primaryVtcs, 0));
        /*
           cout<<"chargediso: "<<mySCFPstruct.chargediso<<endl;
           cout<<"chargediso_primvtx: "<<mySCFPstruct.chargediso_primvtx<<endl;
           cout<<"neutraliso: "<<mySCFPstruct.neutraliso<<endl;
           cout<<"photoniso: "<<mySCFPstruct.photoniso<<endl;
           cout<<"chargediso_rcone: "<<mySCFPstruct.chargediso_rcone<<endl;
           cout<<"chargediso_primvtx_rcone: "<<mySCFPstruct.chargediso_primvtx_rcone<<endl;
           cout<<"neutraliso_rcone: "<<mySCFPstruct.neutraliso_rcone<<endl;
           cout<<"photoniso_rcone: "<<mySCFPstruct.photoniso_rcone<<endl;
           cout<<"eta_rcone: "<<mySCFPstruct.eta_rcone<<endl;
           cout<<"phi_rcone: "<<mySCFPstruct.phi_rcone<<endl;
           cout<<"rcone_isOK: "<<mySCFPstruct.rcone_isOK<<endl;
         */
        TCPhoton::FootprintRemoval foot;
        foot.chargediso = mySCFPstruct.chargediso;
        foot.chargediso_primvtx = mySCFPstruct.chargediso_primvtx;
        foot.neutraliso = mySCFPstruct.neutraliso;
        foot.photoniso = mySCFPstruct.photoniso;
        foot.chargediso_rcone = mySCFPstruct.chargediso_rcone;
        foot.chargediso_primvtx_rcone = mySCFPstruct.chargediso_primvtx_rcone;
        foot.neutraliso_rcone = mySCFPstruct.neutraliso_rcone;
        foot.photoniso_rcone = mySCFPstruct.photoniso_rcone;
        foot.eta_rcone = mySCFPstruct.eta_rcone;
        foot.phi_rcone = mySCFPstruct.phi_rcone;
        foot.rcone_isOK = mySCFPstruct.rcone_isOK;

        myPhoton->SetSCFootprintRemovalStruct(foot);

        //Conversion info
        bool passElectronVeto = !(ConversionTools::hasMatchedPromptElectron(iPhoton->superCluster(), hElectrons, hConversions, vertexBeamSpot.position()));
        myPhoton->SetConversionVeto(passElectronVeto);

        //Effective energy shit
        float phoESEffSigmaRR_x = 0.;
        float phoESEffSigmaRR_y = 0.;
        float phoESEffSigmaRR_z = 0.;

        if (ESRecHits.isValid() && (fabs(iPhoton->superCluster()->eta()) > 1.6 && fabs(iPhoton->superCluster()->eta()) < 3)) {

            vector<float> phoESHits0 = getESHits((*iPhoton).superCluster()->x(), (*iPhoton).superCluster()->y(), (*iPhoton).superCluster()->z(), rechits_map_, geometry_p, topology_p.get(), 0);

            vector<float> phoESShape = getESEffSigmaRR(phoESHits0);
            phoESEffSigmaRR_x = phoESShape[0];
            phoESEffSigmaRR_y = phoESShape[1];
            phoESEffSigmaRR_z = phoESShape[2];
        }
        myPhoton->SetESEffSigmaRR(phoESEffSigmaRR_x, phoESEffSigmaRR_y, phoESEffSigmaRR_z);

        ++photonCount;

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
            }
        }

        //////////////////////
        // Get genParticles //
        //////////////////////

        Handle<GenParticleCollection> genParticleColl;
        iEvent.getByLabel("genParticles", genParticleColl);

        map<const reco::GenParticle*, TCGenParticle*> genMap;
        for (GenParticleCollection::const_iterator myParticle= genParticleColl->begin(); myParticle != genParticleColl->end(); ++myParticle) {

            ////  Leptons and photons and b's, (oh my)
            //// Z's, W's, H's, and now big juicy Gravitons
            if (
                    (abs(myParticle->pdgId()) >= 11 && abs(myParticle->pdgId()) <= 16)
                    || myParticle->pdgId() == 22
                    || abs(myParticle->pdgId()) == 5
                    || abs(myParticle->pdgId()) == 6
                    || abs(myParticle->pdgId()) == 23
                    || abs(myParticle->pdgId()) == 24
                    || abs(myParticle->pdgId()) == 25   // higgs
                    || abs(myParticle->pdgId()) == 35   // another higgs
                    || abs(myParticle->pdgId()) == 36   // more higgses
                    || abs(myParticle->pdgId()) == 39   // graviton (sometimes higgs too)
                    || abs(myParticle->pdgId()) == 443  // jpsi
                    || abs(myParticle->pdgId()) == 553  // upsilon
               ) {
                addGenParticle(&(*myParticle), genPartCount, genMap);

            }
        }



        /////////////////
        // Get genJets //
        /////////////////


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


    ///////////////////
    // Noise filters //
    ///////////////////

    myNoiseFilters.isScraping = false; //isFilteredOutScraping(iEvent, iSetup, 10, 0.25);

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

    Handle<bool> trackingFailureHandle;
    iEvent.getByLabel(trackingFailureTag_, trackingFailureHandle);
    if(trackingFailureHandle.isValid()) myNoiseFilters.isNoiseTracking = !(Bool_t)(*trackingFailureHandle);
    else LogWarning("Filters")<<"tracking Failure NOT valid  ";

    Handle<bool> eeBadScFilterHandle;
    iEvent.getByLabel(eeBadScFilterTag_, eeBadScFilterHandle);
    if(eeBadScFilterHandle.isValid()) myNoiseFilters.isNoiseEEBadSc = !(Bool_t)(*eeBadScFilterHandle);
    else LogWarning("Filters")<<"eeBadSc NOT valid  ";


    Handle<bool> trkPOGFiltersHandle1;
    iEvent.getByLabel(trkPOGFiltersTag1_, trkPOGFiltersHandle1);
    if(trkPOGFiltersHandle1.isValid()) myNoiseFilters.isNoisetrkPOG1 = !(Bool_t)(*trkPOGFiltersHandle1);
    else LogWarning("Filters")<<"trkPOG1 NOT valid  ";

    Handle<bool> trkPOGFiltersHandle2;
    iEvent.getByLabel(trkPOGFiltersTag2_, trkPOGFiltersHandle2);
    if(trkPOGFiltersHandle2.isValid()) myNoiseFilters.isNoisetrkPOG2 = !(Bool_t)(*trkPOGFiltersHandle2);
    else LogWarning("Filters")<<"trkPOG2 NOT valid  ";

    Handle<bool> trkPOGFiltersHandle3;
    iEvent.getByLabel(trkPOGFiltersTag3_, trkPOGFiltersHandle3);
    if(trkPOGFiltersHandle3.isValid()) myNoiseFilters.isNoisetrkPOG3 = !(Bool_t)(*trkPOGFiltersHandle3);
    else LogWarning("Filters")<<"trkPOG3 NOT valid  ";


    edm::Handle<BeamHaloSummary> TheBeamHaloSummary;
    iEvent.getByLabel("BeamHaloSummary",TheBeamHaloSummary);
    const BeamHaloSummary TheSummary = (*TheBeamHaloSummary.product() );
    if (!TheBeamHaloSummary.isValid()) LogWarning("Filters")<<"The Summary (for CSC halo) NOT valid  ";

    myNoiseFilters.isCSCTightHalo = TheSummary.CSCTightHaloId();
    myNoiseFilters.isCSCLooseHalo = TheSummary.CSCLooseHaloId();


    //LogWarning("Filters")<<"\n csc1  "<< myNoiseFilters.isCSCTightHalo<<"  csc2  "<<myNoiseFilters.isCSCLooseHalo
    //  <<" isNoiseHcal HBHE "<<myNoiseFilters.isNoiseHcalHBHE<<"  laser "<<myNoiseFilters.isNoiseHcalLaser<<"\n"
    //  <<" ecal TP  "<<myNoiseFilters.isNoiseEcalTP<<"   ecal BE  "<<myNoiseFilters.isNoiseEcalBE;

    ++nEvents;
    eventTree -> Fill();

    beamSpot->Clear();
    primaryVtx    -> Clear("C");
    recoJets      -> Clear("C");
    recoMuons     -> Clear("C");
    recoElectrons -> Clear("C");
    recoPhotons   -> Clear("C");
    genJets       -> Clear("C");
    genParticles  -> Clear("C");
    triggerObjects.clear();
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
    genJets        = new TClonesArray("TCGenJet");
    genParticles   = new TClonesArray("TCGenParticle");
    beamSpot       = new TVector3();
    recoMET.reset(  new TCMET);

    h1_numOfEvents = fs->make<TH1F>("numOfEvents", "total number of events, unskimmed", 1,0,1);

    eventTree->Branch("recoJets",     &recoJets,       6400, 0);
    eventTree->Branch("recoElectrons",&recoElectrons,  6400, 0);
    eventTree->Branch("recoMuons",    &recoMuons,      6400, 0);
    eventTree->Branch("recoPhotons",  &recoPhotons,    6400, 0);
    eventTree->Branch("recoMET",      recoMET.get(),   6400, 0);
    eventTree->Branch("genJets",      &genJets,        6400, 0);
    eventTree->Branch("genParticles", &genParticles,   6400, 0);

    eventTree->Branch("primaryVtx",      &primaryVtx, 6400, 0);
    eventTree->Branch("beamSpot",        &beamSpot,   6400, 0);
    eventTree->Branch("nPUVertices",     &nPUVertices, "nPUVertices/I");
    eventTree->Branch("nPUVerticesTrue", &nPUVerticesTrue, "nPUVerticesTrue/F");

    eventTree->Branch("isRealData", &isRealData,  "isRealData/O");
    eventTree->Branch("runNumber",  &runNumber,   "runNumber/i");
    eventTree->Branch("eventNumber",&eventNumber, "eventNumber/l");
    eventTree->Branch("lumiSection",&lumiSection, "lumiSection/i");
    eventTree->Branch("bunchCross", &bunchCross,  "bunchCross/i");

    eventTree->Branch("ptHat",      &ptHat,       "ptHat/F");
    eventTree->Branch("qScale",     &qScale,      "qScale/F");
    eventTree->Branch("evtWeight",  &evtWeight,   "evtWeight/F");
    eventTree->Branch("rhoFactor",  &rhoFactor,   "rhoFactor/F");
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

    // Photon PF Iso maker init
    phoIsolator.initializePhotonIsolation(kTRUE);
    phoIsolator.setConeSize(0.3);

    // Electron PF Iso maker init 

    eleIsolator.initializeElectronIsolation(kTRUE);
    eleIsolator.setConeSize(0.4);

    // Initialize Jet PU ID
    myPUJetID.reset(new PileupJetIdAlgo(jetPUIdAlgo_));

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
    //if (isRealData) {
    if (false) {
        edm::Handle<LumiSummary> lumiSummary;
        iLumi.getByLabel("lumiProducer", lumiSummary);

        deliveredLumi  += lumiSummary->avgInsDelLumi()*93.244;
        recordedLumi   += deliveredLumi*lumiSummary->liveFrac();
    }
}


void ntupleProducer::endRun(const edm::Run& iRun, const edm::EventSetup& iSetup)
{
}


void ntupleProducer::endJob()
{
    cout<<nEvents<<endl;
    h1_numOfEvents->SetBinContent(1,nEvents);
    jobTree->Fill();
}


bool ntupleProducer::triggerDecision(edm::Handle<edm::TriggerResults> &hltResults, int iTrigger)
{
    bool triggerPassed = false;
    if(hltResults->wasrun(iTrigger) &&
            hltResults->accept(iTrigger) &&
            !hltResults->error(iTrigger) ){
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

void ntupleProducer::analyzeTrigger(edm::Handle<edm::TriggerResults> &hltResults, edm::Handle<trigger::TriggerEvent> &hltEvent, const std::string& triggerName, int* trigCount) {

    using namespace trigger;

    const unsigned int n(hltConfig_.size());
    const unsigned int triggerIndex(hltConfig_.triggerIndex(triggerName));

    // abort on invalid trigger name
    if(verboseTrigs){
        std::cout<<" n = "<<n<<" triggerIndex = "<<triggerIndex<<" size = "<<hltConfig_.size()<<std::endl;
        std::cout<<" Analyze triggerName : "<<triggerName<<std::endl;
    }

    if (triggerIndex>=n) {
        if(verboseTrigs){
            cout << "DimuonAna::analyzeTrigger: path " << triggerName << " - not found!" << endl;
        }
        return;
    }

    // modules on this trigger path
    const unsigned int moduleIndex(hltResults->index(triggerIndex));
    const unsigned int m(hltConfig_.size(triggerIndex));
    const vector<string>& moduleLabels(hltConfig_.moduleLabels(triggerIndex));

    if (moduleIndex != m-1) return;

    if(verboseTrigs){
        cout << "DimuonAna::analyzeTrigger: path "
            << triggerName << " [" << triggerIndex << "]" << endl;

        std::cout<<" n = "<< n<<" triggerIndex = "<<triggerIndex<<" m = "<<m<<std::endl;
        std::cout<<" moduleLabels = "<<moduleLabels.size()<<" moduleIndex = "<<moduleIndex<<std::endl;

        // Results from TriggerResults product
        cout << " Trigger path status:"
            << " WasRun=" << hltResults->wasrun(triggerIndex)
            << " Accept=" << hltResults->accept(triggerIndex)
            << " Error =" << hltResults->error(triggerIndex)
            << endl;
        cout << " Last active module - label/type: "
            << moduleLabels[moduleIndex] << "/" << hltConfig_.moduleType(moduleLabels[moduleIndex])
            << " [" << moduleIndex << " out of 0-" << (m-1) << " on this path]"
            << endl;
    }
    assert(moduleIndex<m);

    // Results from TriggerEvent product - Attention: must look only for
    // modules actually run in this path for this event!
    //std::vector < GlobalVector > passMomenta;
    for (unsigned int j=0; j<=moduleIndex; ++j) {
        const string&   moduleLabel(moduleLabels[j]);
        const string    moduleType(hltConfig_.moduleType(moduleLabel));

        // check whether the module is packed up in TriggerEvent product
        //cout<<hltEvent->filterIndex(InputTag(moduleLabel,"",hltProcess_))<<endl;

        const unsigned int filterIndex(hltEvent->filterIndex(InputTag(moduleLabel,"",hltProcess_)));

        // if ( (moduleLabel.find("Calo") == string::npos) )continue;
        // if ( (moduleLabel.find("hltEventle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoDZ") == string::npos)
        // && (moduleLabel.find("hltEventle17CaloId") == string::npos)
        // && (moduleLabel.find("hltEventle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoDoubleFilter") == string::npos) ) continue;

        if(verboseTrigs){
            std::cout<<" j = "<<j<<" modLabel/moduleType = "<<moduleLabel<<"/"<<moduleType<<" filterIndex = "<<filterIndex<<" sizeF = "<<hltEvent->sizeFilters()<<std::endl;
        }

        if (filterIndex < hltEvent->sizeFilters()) {
            if(verboseTrigs){
                cout << " 'L3' (or 'L1', 'L2') filter in slot " << j << " - label/type " << moduleLabel << "/" << moduleType << endl;
            }
            const Vids& VIDS (hltEvent->filterIds(filterIndex));
            const Keys& KEYS(hltEvent->filterKeys(filterIndex));
            const size_type nI(VIDS.size());
            const size_type nK(KEYS.size());
            assert(nI==nK);
            const size_type n(max(nI,nK));
            if(verboseTrigs){
                cout << " " << n << " accepted 'L3' (or 'L1', 'L2') objects found: " << endl;
            }
            const TriggerObjectCollection& TOC(hltEvent->getObjects());
            for (size_type i = 0; i != n; ++i) {
                const TriggerObject& TO(TOC[KEYS[i]]);
                if (TO.pt() < 5 || (fabs(TO.id()) != 13 && fabs(TO.id()) != 11)) continue;
                TCTriggerObject trigObj = TCTriggerObject();

                //cout<<" i = "<<i<<" moduleLabel/moduleType : "<<moduleLabel<<"/"<<moduleType<<endl;
                //cout << " " << i << " " << VIDS[i] << "/" << KEYS[i] << ": \n"
                // << "triggerName = "<<triggerName<<" moduleLabel="<<moduleLabel<<"\n "
                // << TO.id() << " " << TO.pt() << " " << TO.eta() << " " << TO.phi() << " " << TO.mass() << " " << TO.energy()
                //<< endl;

                if (TO.mass() <= 0 || TO.pt() <= 0) continue;

                trigObj.SetPtEtaPhiE(TO.pt(), TO.eta(), TO.phi(), TO.energy());
                trigObj.SetHLTName(triggerName);
                trigObj.SetModuleName(moduleLabel);
                trigObj.SetId(TO.id());

                float minDeltaR = 99.;
                for (unsigned j = 0; j < triggerObjects.size(); ++j) {
                    float deltaR = triggerObjects[j].DeltaR(trigObj);
                    if (deltaR < minDeltaR && trigObj.GetId() == triggerObjects[j].GetId()) minDeltaR = deltaR;
                }

                if (minDeltaR < 0.1) continue;

                triggerObjects.push_back(trigObj);
                (*trigCount)++;
            }
        }
    }
    //cout<<endl;
}

float ntupleProducer::MatchBTagsToJets(const reco::JetTagCollection bTags,const reco::PFJet jet)
{
    float discrValue = -999;
    for (tag_iter iTag = bTags.begin(); iTag != bTags.end(); iTag++) {
        if (sqrt(pow(iTag->first->eta() - jet.eta(), 2) + pow(deltaPhi(iTag->first->phi(),jet.phi()), 2)) == 0.) {
            discrValue = iTag->second;
            break;
        }
    }

    return discrValue;
}

void ntupleProducer::initJetEnergyCorrector(const edm::EventSetup &iSetup, bool isData)
{
    //jet energy correction levels to apply on raw jet
    std::vector<JetCorrectorParameters> jetCorPars_;
    std::vector<std::string> jecLevels;
    jecLevels.push_back("L1FastJet");
    jecLevels.push_back("L2Relative");
    jecLevels.push_back("L3Absolute");
    if(isData) jecLevels.push_back("L2L3Residual");

    //check the corrector parameters needed according to the correction levels
    edm::ESHandle<JetCorrectorParametersCollection> parameters;
    iSetup.get<JetCorrectionsRecord>().get(jecTag_,parameters);
    for(std::vector<std::string>::const_iterator ll = jecLevels.begin(); ll != jecLevels.end(); ++ll)
    {
        const JetCorrectorParameters& ip = (*parameters)[*ll];
        jetCorPars_.push_back(ip);
    }

    //instantiate the jet corrector
    jecCor.reset(new FactorizedJetCorrector(jetCorPars_));
}


TCGenParticle* ntupleProducer::addGenParticle(const reco::GenParticle* myParticle, int& genPartCount, map<const reco::GenParticle*,TCGenParticle*>& genMap)
{
    TCGenParticle* genCon;
    map<const reco::GenParticle*,TCGenParticle*>::iterator it;
    it = genMap.find(myParticle);
    if (it == genMap.end()){
        genCon = new ((*genParticles)[genPartCount]) TCGenParticle;
        ++genPartCount;
        genMap[myParticle] = genCon;
        genCon->SetPxPyPzE(myParticle->px(), myParticle->py(), myParticle->pz(), myParticle->energy() );
        genCon->SetVtx(myParticle->vx(), myParticle->vy(), myParticle->vz());
        genCon->SetCharge(myParticle->charge());
        genCon->SetPDGId( myParticle->pdgId());
        genCon->SetStatus(myParticle->status());
        map<const reco::GenParticle*,TCGenParticle*>::iterator momIt;

        genCon->SetMotherId(myParticle->mother()->pdgId());

        if (myParticle->numberOfMothers() == 0){
            genCon->SetMother(0);
        }else if(
                abs(myParticle->mother()->pdgId()) != 6
                && abs(myParticle->mother()->pdgId()) != 5
                && abs(myParticle->mother()->pdgId()) != 11
                && abs(myParticle->mother()->pdgId()) != 12
                && abs(myParticle->mother()->pdgId()) != 13
                && abs(myParticle->mother()->pdgId()) != 14
                && abs(myParticle->mother()->pdgId()) != 15
                && abs(myParticle->mother()->pdgId()) != 16
                && abs(myParticle->mother()->pdgId()) != 22
                && abs(myParticle->mother()->pdgId()) != 23
                && abs(myParticle->mother()->pdgId()) != 24
                && abs(myParticle->mother()->pdgId()) != 25
                && abs(myParticle->mother()->pdgId()) != 35
                && abs(myParticle->mother()->pdgId()) != 36
                && abs(myParticle->mother()->pdgId()) != 39
                && abs(myParticle->mother()->pdgId()) != 443  //Jpsi
                && abs(myParticle->mother()->pdgId()) != 553  //Upsilon
                )
        {
            genCon->SetMother(0);
        }else{
            momIt = genMap.find((const reco::GenParticle*)myParticle->mother());
            if (momIt == genMap.end()){
                genCon->SetMother(addGenParticle((const reco::GenParticle*)myParticle->mother(), genPartCount, genMap));
            }else{
                genCon->SetMother(momIt->second);
            }
        }
    }
    else
        genCon = it->second;

    return genCon;
}

vector<float> ntupleProducer::getESHits(double X, double Y, double Z, map<DetId, EcalRecHit> rechits_map, const CaloSubdetectorGeometry*& geometry_p, CaloSubdetectorTopology *topology_p, int row) {

    //cout<<row<<endl;

    vector<float> esHits;

    //double X = bcPtr->x();
    //double Y = bcPtr->y();
    //double Z = bcPtr->z();
    const GlobalPoint point(X,Y,Z);

    DetId esId1 = (dynamic_cast<const EcalPreshowerGeometry*>(geometry_p))->getClosestCellInPlane(point, 1);
    DetId esId2 = (dynamic_cast<const EcalPreshowerGeometry*>(geometry_p))->getClosestCellInPlane(point, 2);
    ESDetId esDetId1 = (esId1 == DetId(0)) ? ESDetId(0) : ESDetId(esId1);
    ESDetId esDetId2 = (esId2 == DetId(0)) ? ESDetId(0) : ESDetId(esId2);

    map<DetId, EcalRecHit>::iterator it;
    ESDetId next;
    ESDetId strip1;
    ESDetId strip2;

    strip1 = esDetId1;
    strip2 = esDetId2;

    EcalPreshowerNavigator theESNav1(strip1, topology_p);
    theESNav1.setHome(strip1);

    EcalPreshowerNavigator theESNav2(strip2, topology_p);
    theESNav2.setHome(strip2);

    if (row == 1) {
        if (strip1 != ESDetId(0)) strip1 = theESNav1.north();
        if (strip2 != ESDetId(0)) strip2 = theESNav2.east();
    } else if (row == -1) {
        if (strip1 != ESDetId(0)) strip1 = theESNav1.south();
        if (strip2 != ESDetId(0)) strip2 = theESNav2.west();
    }

    // Plane 2
    if (strip1 == ESDetId(0)) {
        for (unsigned int i=0; i<31; ++i) esHits.push_back(0);
    } else {

        it = rechits_map.find(strip1);
        if (it->second.energy() > 1.0e-10 && it != rechits_map.end()) esHits.push_back(it->second.energy());
        else esHits.push_back(0);
        //cout<<"center : "<<strip1<<" "<<it->second.energy()<<endl;

        // east road
        for (unsigned int i=0; i<15; ++i) {
            next = theESNav1.east();
            if (next != ESDetId(0)) {
                it = rechits_map.find(next);
                if (it->second.energy() > 1.0e-10 && it != rechits_map.end()) esHits.push_back(it->second.energy());
                else esHits.push_back(0);
                //cout<<"east "<<i<<" : "<<next<<" "<<it->second.energy()<<endl;
            } else {
                for (unsigned int j=i; j<15; ++j) esHits.push_back(0);
                break;
                //cout<<"east "<<i<<" : "<<next<<" "<<0<<endl;
            }
        }

        // west road
        theESNav1.setHome(strip1);
        theESNav1.home();
        for (unsigned int i=0; i<15; ++i) {
            next = theESNav1.west();
            if (next != ESDetId(0)) {
                it = rechits_map.find(next);
                if (it->second.energy() > 1.0e-10 && it != rechits_map.end()) esHits.push_back(it->second.energy());
                else esHits.push_back(0);
                //cout<<"west "<<i<<" : "<<next<<" "<<it->second.energy()<<endl;
            } else {
                for (unsigned int j=i; j<15; ++j) esHits.push_back(0);
                break;
                //cout<<"west "<<i<<" : "<<next<<" "<<0<<endl;
            }
        }
    }

    if (strip2 == ESDetId(0)) {
        for (unsigned int i=0; i<31; ++i) esHits.push_back(0);
    } else {

        it = rechits_map.find(strip2);
        if (it->second.energy() > 1.0e-10 && it != rechits_map.end()) esHits.push_back(it->second.energy());
        else esHits.push_back(0);
        //cout<<"center : "<<strip2<<" "<<it->second.energy()<<endl;

        // north road
        for (unsigned int i=0; i<15; ++i) {
            next = theESNav2.north();
            if (next != ESDetId(0)) {
                it = rechits_map.find(next);
                if (it->second.energy() > 1.0e-10 && it != rechits_map.end()) esHits.push_back(it->second.energy());
                else esHits.push_back(0);
                //cout<<"north "<<i<<" : "<<next<<" "<<it->second.energy()<<endl;
            } else {
                for (unsigned int j=i; j<15; ++j) esHits.push_back(0);
                break;
                //cout<<"north "<<i<<" : "<<next<<" "<<0<<endl;
            }
        }

        // south road
        theESNav2.setHome(strip2);
        theESNav2.home();
        for (unsigned int i=0; i<15; ++i) {
            next = theESNav2.south();
            if (next != ESDetId(0)) {
                it = rechits_map.find(next);
                if (it->second.energy() > 1.0e-10 && it != rechits_map.end()) esHits.push_back(it->second.energy());
                else esHits.push_back(0);
                //cout<<"south "<<i<<" : "<<next<<" "<<it->second.energy()<<endl;
            } else {
                for (unsigned int j=i; j<15; ++j) esHits.push_back(0);
                break;
                //cout<<"south "<<i<<" : "<<next<<" "<<0<<endl;
            }
        }
    }

    return esHits;
}

vector<float> ntupleProducer::getESEffSigmaRR(vector<float> ESHits0)
{
    const int nBIN = 21;
    vector<float> esShape;

    TH1F *htmpF = new TH1F("htmpF","",nBIN,0,nBIN);
    TH1F *htmpR = new TH1F("htmpR","",nBIN,0,nBIN);
    htmpF->Reset(); htmpR->Reset();

    Float_t effsigmaRR=0.;

    for(int ibin=0; ibin<((nBIN+1)/2); ++ibin) {
        if (ibin==0) {
            htmpF->SetBinContent((nBIN+1)/2,ESHits0[ibin]);
            htmpR->SetBinContent((nBIN+1)/2,ESHits0[ibin+31]);
        } else { // hits sourd the seed
            htmpF->SetBinContent((nBIN+1)/2+ibin,ESHits0[ibin]);
            htmpF->SetBinContent((nBIN+1)/2-ibin,ESHits0[ibin+15]);
            htmpR->SetBinContent((nBIN+1)/2+ibin,ESHits0[ibin+31]);
            htmpR->SetBinContent((nBIN+1)/2-ibin,ESHits0[ibin+31+15]);
        }
    }

    // ---- Effective Energy Deposit Width ---- //
    double EffWidthSigmaXX = 0.;
    double EffWidthSigmaYY = 0.;
    double totalEnergyXX   = 0.;
    double totalEnergyYY   = 0.;
    double EffStatsXX      = 0.;
    double EffStatsYY      = 0.;
    for (int id_X=1; id_X<=21; ++id_X) {
        totalEnergyXX  += htmpF->GetBinContent(id_X);
        EffStatsXX     += htmpF->GetBinContent(id_X)*(id_X-11)*(id_X-11);
        totalEnergyYY  += htmpR->GetBinContent(id_X);
        EffStatsYY     += htmpR->GetBinContent(id_X)*(id_X-11)*(id_X-11);
    }
    // If denominator == 0, effsigmaRR = 0;
    EffWidthSigmaXX  = (totalEnergyXX>0.)  ? sqrt(fabs(EffStatsXX  / totalEnergyXX))   : 0.;
    EffWidthSigmaYY  = (totalEnergyYY>0.)  ? sqrt(fabs(EffStatsYY  / totalEnergyYY))   : 0.;
    effsigmaRR =  ((totalEnergyXX  + totalEnergyYY) >0.) ? sqrt(EffWidthSigmaXX  * EffWidthSigmaXX  + EffWidthSigmaYY  * EffWidthSigmaYY)  : 0.;
    esShape.push_back(effsigmaRR);
    esShape.push_back(EffWidthSigmaXX);
    esShape.push_back(EffWidthSigmaYY);

    delete htmpF;
    delete htmpR;

    return esShape;
}


TCTrack::ConversionInfo ntupleProducer::CheckForConversions(const edm::Handle<reco::ConversionCollection> &convCol,
        const reco::GsfTrackRef &gsf,
        const math::XYZPoint &bs, const math::XYZPoint &pv)
{
    TCTrack::ConversionInfo * convInfo = new TCTrack::ConversionInfo();
    //int iconv=-1;
    for (reco::ConversionCollection::const_iterator conv = convCol->begin(); conv!= convCol->end(); ++conv) {
        //iconv++;

        reco::Vertex vtx = conv->conversionVertex();
        if (vtx.isValid()) {
            if (ConversionTools::matchesConversion(gsf, *conv)) {

                (*convInfo).isValid = true;

                (*convInfo).vtxProb = TMath::Prob( vtx.chi2(), vtx.ndof() );
                math::XYZVector mom(conv->refittedPairMomentum());
                double dbsx = vtx.x() - bs.x();
                double dbsy = vtx.y() - bs.y();
                (*convInfo).lxyBS = (mom.x()*dbsx + mom.y()*dbsy)/mom.rho();

                double dpvx = vtx.x() - pv.x();
                double dpvy = vtx.y() - pv.y();
                (*convInfo).lxyPV = (mom.x()*dpvx + mom.y()*dpvy)/mom.rho();

                (*convInfo).nHitsMax=0;
                for (std::vector<uint8_t>::const_iterator it = conv->nHitsBeforeVtx().begin(); it!=conv->nHitsBeforeVtx().end(); ++it) {
                    if ((*it)>(*convInfo).nHitsMax) (*convInfo).nHitsMax = (*it);
                }

                break;
            }
        }
    }
    return (*convInfo);
}

void ntupleProducer::MatchTriggerObject(TCPhysObject& physObj, const unsigned pdgID)
{
    for (unsigned j = 0; j < triggerObjects.size(); ++j) {
        float deltaR    = triggerObjects[j].DeltaR(physObj);
        //float deltaPt   = fabs(triggerObjects[j].Pt() - physObj.Pt())/physObj.Pt();
        if (deltaR < 0.1 && fabs(triggerObjects[j].GetId()) == pdgID) {
            //cout << triggerObjects[j].GetHLTName() << "\t" << triggerObjects[j].GetModuleName() << "\t" << deltaR << "\t" << deltaPt << endl;;
            physObj.SetTriggered(true);
            break;
        } 
    }
}

//define this as a plug-in
DEFINE_FWK_MODULE(ntupleProducer);

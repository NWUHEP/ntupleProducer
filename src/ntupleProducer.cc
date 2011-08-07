#include "ntupleProducer.h"

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

	int vtxCount, jetCount, metCount, muCount, eleCount, photonCount, tauCount, genCount, genPartCount;
	vtxCount = jetCount = metCount = muCount = eleCount = photonCount = tauCount = genCount = genPartCount = 0;
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

		edm::Handle<reco::JetTagCollection> bTagCollectionTCHE;
		iEvent.getByLabel("trackCountingHighEffBJetTags", bTagCollectionTCHE);
		const reco::JetTagCollection & bTagsTCHE = *(bTagCollectionTCHE.product());

		edm::Handle<reco::JetTagCollection> bTagCollectionSSV;
		iEvent.getByLabel("simpleSecondaryVertexHighEffBJetTags", bTagCollectionSSV);
		const reco::JetTagCollection & bTagsSSV = *(bTagCollectionSSV.product());
		typedef reco::JetTagCollection::const_iterator tag_iter;

		edm::ESHandle<JetCorrectorParametersCollection> JetCorParColl;
		iSetup.get<JetCorrectionsRecord>().get("AK5PF",JetCorParColl);
		JetCorrectorParameters const & JetCorPar = (*JetCorParColl)["Uncertainty"];
		JetCorrectionUncertainty *jecUncertainty = new JetCorrectionUncertainty(JetCorPar);

		const JetCorrector* correctorL1  = JetCorrector::getJetCorrector("ak5PFL1Fastjet",iSetup);
		const JetCorrector* correctorL2  = JetCorrector::getJetCorrector("ak5PFL2Relative",iSetup);
		const JetCorrector* correctorL3  = JetCorrector::getJetCorrector("ak5PFL3Absolute",iSetup);
		const JetCorrector* correctorRes = JetCorrector::getJetCorrector("ak5PFResidual", iSetup);

		Handle<reco::PFJetCollection> PFJets;
		iEvent.getByLabel(jetTag_, PFJets);

		for (PFJetCollection::const_iterator jet_iter = PFJets->begin(); jet_iter!= PFJets->end(); ++jet_iter) {

			reco::PFJet myJet  = reco::PFJet(*jet_iter);
			reco::PFJet corJet = reco::PFJet(*jet_iter);

			int index = jet_iter - PFJets->begin();
			edm::RefToBase<reco::Jet> jetRef(edm::Ref<reco::PFJetCollection>(PFJets,index));

			float scale1 = correctorL1->correction(corJet, jetRef, iEvent, iSetup);
			corJet.scaleEnergy(scale1);
			float scale2 = correctorL2->correction(corJet);
			corJet.scaleEnergy(scale2);
			float scale3 = correctorL3->correction(corJet);
			corJet.scaleEnergy(scale3);

			if (corJet.pt() < 10.) continue;

			TCJet* jetCon = new ((*recoJets)[jetCount]) TCJet;

			jetCon->SetP4(myJet.px(), myJet.py(), myJet.pz(), myJet.energy());
			jetCon->SetVtx(0., 0., 0.);
			jetCon->SetChHadFrac(myJet.chargedHadronEnergyFraction());
			jetCon->SetNeuHadFrac(myJet.neutralHadronEnergyFraction());
			jetCon->SetChEmFrac(myJet.chargedEmEnergyFraction());
			jetCon->SetNeuEmFrac(myJet.neutralEmEnergyFraction());
			jetCon->SetNumConstit(myJet.chargedMultiplicity() + myJet.neutralMultiplicity());
			jetCon->SetNumChPart(myJet.chargedMultiplicity());

			// Get b-tag information
			for (tag_iter iTag = bTagsTCHE.begin(); iTag != bTagsTCHE.end(); iTag++) {
				if (sqrt(pow(iTag->first->eta() - corJet.eta(), 2) + pow(deltaPhi(iTag->first->phi(),corJet.phi()), 2)) == 0.) {
					jetCon->SetBDiscrTrkCountHiEff(iTag->second);
				}
			}
			for (tag_iter iTag = bTagsSSV.begin(); iTag != bTagsSSV.end(); iTag++) {
				if (sqrt(pow(iTag->first->eta() - corJet.eta(), 2) + pow(deltaPhi(iTag->first->phi(),corJet.phi()), 2)) == 0.) {
					jetCon->SetBDiscrSecVtxSimple(iTag->second);
				}
			}

			jetCon->SetJetCorr(1, scale1);
			jetCon->SetJetCorr(2, scale2);
			jetCon->SetJetCorr(3, scale3);

			if (isRealData) {
				float scaleRes = correctorRes->correction(corJet);
				jetCon->SetJetCorr(4, scaleRes);
			}

			jecUncertainty->setJetEta(corJet.eta());
			jecUncertainty->setJetPt(corJet.pt());
			jetCon->SetUncertaintyJES(jecUncertainty->getUncertainty(true)); 

			/////////////////////////
			// Associate to vertex //
			/////////////////////////

			if(fabs(myJet.eta()) < 2.5){
				associateJetToVertex(myJet, primaryVtcs, jetCon);
			} else {
				jetCon->SetVtxSumPtFrac(-1);
				jetCon->SetVtxSumPt(-1);
				jetCon->SetVtxTrackFrac(-1);
				jetCon->SetVtxNTracks(-1);
				jetCon->SetVtxIndex(0);
			}
			++jetCount;
		}   
		delete jecUncertainty;
	}



	/////////////
	// Get MET //
	/////////////

	if (saveMET_){

		Handle<PFMETCollection> MET;
		iEvent.getByLabel(metTag_, MET);
		PFMETCollection::const_iterator pfMET = MET->begin();

		if (MET.isValid()) {
			recoMET->SetSumEt(pfMET->sumEt());
			recoMET->SetMet(pfMET->et());
			recoMET->SetPhi(pfMET->phi());
			recoMET->SetPhotonEtFraction(pfMET->photonEtFraction());
			recoMET->SetElectronEtFraction(pfMET->electronEtFraction());
			recoMET->SetMuonEtFraction(pfMET->muonEtFraction());
			recoMET->SetNeutralHadronEtFraction(pfMET->neutralHadronEtFraction());
			recoMET->SetChargedHadronEtFraction(pfMET->chargedHadronEtFraction());
			recoMET->SetHFHadronEtFraction(pfMET->HFHadronEtFraction());
			recoMET->SetHFEMEtFraction(pfMET->HFEMEtFraction());

			Handle<PFMETCollection> corMET;
			iEvent.getByLabel("metJESCorAK5PF", corMET);
			reco::PFMET iMET = corMET->front();
			recoMET->SetCorrectedSumEt(iMET.sumEt());
			recoMET->SetCorrectedMet(iMET.et());
			recoMET->SetCorrectedPhi(iMET.phi());
		}
	}

	///////////////
	// Get muons //
	///////////////

	if (saveMuons_) {

		Handle<MuonCollection> muons;
		iEvent.getByLabel(muonTag_, muons);

		for (MuonCollection::const_iterator mu = muons->begin(); mu != muons->end(); ++mu) {
			if (!(mu->isGlobalMuon() && mu->isTrackerMuon()) || mu->pt() < 10.) continue;
			TCMuon* muCon = new ((*recoMuons)[muCount]) TCMuon;

			muCon->SetP4(mu->px(), mu->py(), mu->pz(), mu->energy());
			muCon->SetVtx(mu->globalTrack()->vx(),mu->globalTrack()->vy(),mu->globalTrack()->vz());
			muCon->SetPtError(mu->globalTrack()->ptError());
			muCon->SetCharge(mu->charge());
			muCon->SetIsGLB(mu->isGlobalMuon());
			muCon->SetIsTRK(mu->isTrackerMuon());
			muCon->SetNumberOfMatches(mu->numberOfMatches());
			muCon->SetNumberOfValidPixelHits(mu->globalTrack()->hitPattern().numberOfValidPixelHits());
			muCon->SetNumberOfValidTrackerHits(mu->globalTrack()->hitPattern().numberOfValidTrackerHits()); 
			muCon->SetNumberOfValidMuonHits(mu->globalTrack()->hitPattern().numberOfValidMuonHits());
			muCon->SetNumberOfLostPixelHits(mu->globalTrack()->hitPattern().numberOfLostPixelHits());
			muCon->SetNumberOfLostTrackerHits(mu->globalTrack()->hitPattern().numberOfLostTrackerHits());
			muCon->SetNormalizedChi2(mu->globalTrack()->normalizedChi2());

			muCon->SetCaloComp(mu->caloCompatibility());
			muCon->SetSegComp(muon::segmentCompatibility(*mu));

			muCon->SetNtracks03(mu->isolationR03().nTracks);
			muCon->SetEmIso03(mu->isolationR03().emEt);
			muCon->SetHadIso03(mu->isolationR03().hadEt);
			muCon->SetTrkIso03(mu->isolationR03().sumPt);

			muCon->SetNtracks05(mu->isolationR05().nTracks);
			muCon->SetEmIso05(mu->isolationR05().emEt);
			muCon->SetHadIso05(mu->isolationR05().hadEt);
			muCon->SetTrkIso05(mu->isolationR05().sumPt);

			float sumPt3 = 0;
			float gamma3 = 0;
			float neutral3 = 0;
			float sumPt4 = 0;
			float gamma4 = 0;
			float neutral4 = 0;
			float sumPt5 = 0;
			float gamma5 = 0;
			float neutral5 = 0;

			muCon->SetPfSumPt(0.5, sumPt5);
			muCon->SetPfSumPt(0.4, sumPt4);
			muCon->SetPfSumPt(0.3, sumPt3);
			muCon->SetPfEGamma(0.5, gamma5);
			muCon->SetPfEGamma(0.4, gamma4);
			muCon->SetPfEGamma(0.3, gamma3);
			muCon->SetPfENeutral(0.5, neutral5);
			muCon->SetPfENeutral(0.4, neutral4);
			muCon->SetPfENeutral(0.3, neutral3);

			muCount++;
		}
	}


	///////////////////
	// Get electrons //
	///////////////////

	if (saveElectrons_) {

		Handle<GsfElectronCollection> electrons;
		iEvent.getByLabel(electronTag_, electrons);

		edm::Handle<edm::ValueMap<float> > eIDValueMap95;
		iEvent.getByLabel( "simpleEleId95relIso" , eIDValueMap95 );
		const edm::ValueMap<float> & eIDmap95 = * eIDValueMap95 ;

		edm::Handle<edm::ValueMap<float> > eIDValueMap90;
		iEvent.getByLabel( "simpleEleId90relIso" , eIDValueMap90 );
		const edm::ValueMap<float> & eIDmap90 = * eIDValueMap90 ;

		edm::Handle<edm::ValueMap<float> > eIDValueMap85;
		iEvent.getByLabel( "simpleEleId85relIso" , eIDValueMap85 );
		const edm::ValueMap<float> & eIDmap85 = * eIDValueMap85 ;

		edm::Handle<edm::ValueMap<float> > eIDValueMap80;
		iEvent.getByLabel( "simpleEleId80relIso" , eIDValueMap80 );
		const edm::ValueMap<float> & eIDmap80 = * eIDValueMap80 ;

		edm::Handle<edm::ValueMap<float> > eIDValueMap70;
		iEvent.getByLabel( "simpleEleId70relIso" , eIDValueMap70 );
		const edm::ValueMap<float> & eIDmap70 = * eIDValueMap70 ;

		edm::Handle<edm::ValueMap<float> > eIDValueMap60;
		iEvent.getByLabel( "simpleEleId60relIso" , eIDValueMap60 );
		const edm::ValueMap<float> & eIDmap60 = * eIDValueMap60;

		for (unsigned int i = 0; i < electrons->size(); i++) {
			edm::Ref<reco::GsfElectronCollection> electronRef(electrons,i);

			if (electronRef->pt() < 10) continue;

			int cuts95 = eIDmap95[electronRef];
			int cuts90 = eIDmap90[electronRef];
			int cuts85 = eIDmap85[electronRef];
			int cuts80 = eIDmap80[electronRef];
			int cuts70 = eIDmap70[electronRef];
			int cuts60 = eIDmap60[electronRef];

			TCElectron* eleCon = new ((*recoElectrons)[eleCount]) TCElectron;

			eleCon->SetP4(electronRef->px(),electronRef->py(),electronRef->pz(),electronRef->p());
			eleCon->SetVtx(electronRef->gsfTrack()->vx(),electronRef->gsfTrack()->vy(),electronRef->gsfTrack()->vz());
			eleCon->SetCharge(electronRef->charge());

			eleCon->SetNumberOfValidPixelHits(electronRef->gsfTrack()->hitPattern().numberOfValidPixelHits());
			eleCon->SetNumberOfValidTrackerHits(electronRef->gsfTrack()->hitPattern().numberOfValidTrackerHits());
			eleCon->SetNumberOfLostPixelHits(electronRef->gsfTrack()->hitPattern().numberOfLostPixelHits());
			eleCon->SetNumberOfLostTrackerHits(electronRef->gsfTrack()->hitPattern().numberOfLostTrackerHits());

			eleCon->SetIsEB(electronRef->isEB());
			eleCon->SetIsEE(electronRef->isEE());
			eleCon->SetIsInGap(electronRef->isGap());

			eleCon->SetEmIso03( electronRef->dr03EcalRecHitSumEt());
			eleCon->SetHadIso03(electronRef->dr03HcalTowerSumEt());
			eleCon->SetTrkIso03(electronRef->dr03TkSumPt());
			eleCon->SetEmIso04( electronRef->dr04EcalRecHitSumEt());
			eleCon->SetHadIso04(electronRef->dr04HcalTowerSumEt());
			eleCon->SetTrkIso04(electronRef->dr04TkSumPt());

			eleCon->SetHadOverEm(electronRef->hadronicOverEm());
			eleCon->SetDphiSuperCluster(electronRef->deltaPhiSuperClusterTrackAtVtx());
			eleCon->SetDetaSuperCluster(electronRef->deltaEtaSuperClusterTrackAtVtx());
			eleCon->SetSigmaIetaIeta(electronRef->sigmaIetaIeta());

			eleCon->SetConversionFlag(electronRef->convFlags());
			eleCon->SetConversionDist(electronRef->convDist());
			eleCon->SetConversionDcot(electronRef->convDcot());
			eleCon->SetConversionRad(electronRef->convRadius());

			eleCon->SetCutLevel(cuts95, 95);
			eleCon->SetCutLevel(cuts90, 90);
			eleCon->SetCutLevel(cuts85, 85);
			eleCon->SetCutLevel(cuts80, 80);
			eleCon->SetCutLevel(cuts70, 70);
			eleCon->SetCutLevel(cuts60, 60);

			//float sumPt3 = 0;
			//float gamma3 = 0;
			//float neutral3 = 0;
			//float sumPt4 = 0;
			//float gamma4 = 0;
			//float neutral4 = 0;
			//float sumPt5 = 0;
			//float gamma5 = 0;
			//float neutral5 = 0;
			//eleCon->SetPfSumPt(0.5, sumPt5);
			//eleCon->SetPfEGamma(0.5, gamma5);
			//eleCon->SetPfSumPt(0.4, sumPt4);
			//eleCon->SetPfEGamma(0.4, gamma4);
			//eleCon->SetPfENeutral(0.5, neutral5);
			//eleCon->SetPfENeutral(0.4, neutral4);
			eleCon->SetPfEGamma(0.3, electronRef->pfIsolationVariables().photonIso);
			eleCon->SetPfSumPt(0.3, electronRef->pfIsolationVariables().chargedHadronIso);
			eleCon->SetPfENeutral(0.3, electronRef->pfIsolationVariables().neutralHadronIso);

			eleCount++;
		}
	}

	/////////////////
	// Get photons //
	/////////////////

	if (savePhotons_) {
		edm::Handle<reco::PhotonCollection> photons;
		iEvent.getByLabel(photonTag_, photons);

		for (reco::PhotonCollection::const_iterator iPhoton = photons->begin(); iPhoton != photons->end() ; ++iPhoton)
		{
            
			TCPhoton* myPhoton = new ((*recoPhotons)[photonCount]) TCPhoton;
			myPhoton->SetP4(iPhoton->px(), iPhoton->py(), iPhoton->pz(), iPhoton->p());
			myPhoton->SetEMIso(iPhoton->ecalRecHitSumEtConeDR04());
			myPhoton->SetHADIso(iPhoton->hcalTowerSumEtConeDR04());
			myPhoton->SetTRKIso(iPhoton->trkSumPtHollowConeDR04());
			myPhoton->SetHadOverEm(iPhoton->hadronicOverEm());
			myPhoton->SetSigmaIEtaIEta(iPhoton->sigmaIetaIeta());
            myPhoton->SetR9(iPhoton->r9());
			//myPhoton->SetSigmaIphiIphi();
			//myPhoton->SetE2OverE9();
			myPhoton->SetEtaSupercluster(iPhoton->superCluster()->eta());
			myPhoton->SetTrackVeto(iPhoton->hasPixelSeed());

			++photonCount;
		}
	}

	//////////////
	// Get taus //
	//////////////


	////////////////////////
	// Get gen-level info //
	////////////////////////


	if (!isRealData) {

		Handle<GenEventInfoProduct> GenEventInfoHandle;
		iEvent.getByLabel("generator", GenEventInfoHandle);

		evtWeight = ptHat = qScale = -1;

		if (GenEventInfoHandle.isValid()) {
			//qScale       = GenEventInfoHandle->qScale();
			//ptHat        = (GenEventInfoHandle->hasBinningValues() ? GenEventInfoHandle->binningValues()[0] : 0.0);
			//evtWeight    = GenEventInfoHandle->weight();
		}


		////////////////////
		// PU information //
		////////////////////


		Handle<std::vector< PileupSummaryInfo > >  PUInfo;
		iEvent.getByLabel(edm::InputTag("addPileupInfo"), PUInfo);
		std::vector<PileupSummaryInfo>::const_iterator iPV;

		for(iPV = PUInfo->begin(); iPV != PUInfo->end(); ++iPV) if (iPV->getBunchCrossing() == 0) nPUVertices = iPV->getPU_NumInteractions();


		//////////////////////
		// Get genParticles //
		//////////////////////

        Handle<GenParticleCollection> genParticleColl;
        iEvent.getByLabel("genParticles", genParticleColl);

        for (GenParticleCollection::const_iterator iGenPart = genParticleColl->begin(); iGenPart != genParticleColl->end(); ++iGenPart) {
            const reco::GenParticle myParticle = reco::GenParticle(*iGenPart);

            if (myParticle.status() == 3 && (abs(myParticle.pdgId()) == 6 || abs(myParticle.pdgId()) == 23 || abs(myParticle.pdgId()) == 24)) {
                for (size_t i = 0; i < myParticle.numberOfDaughters(); ++i) {
                    const reco::Candidate *myDaughter = myParticle.daughter(i);
                    if (abs(myDaughter->pdgId()) == 5 || (abs(myDaughter->pdgId()) >= 11 && abs(myDaughter->pdgId()) <= 16)) {
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

			edm::Handle<reco::JetFlavourMatchingCollection> jetFlavourMC;
			iEvent.getByLabel("GenJetbyValAlgo", jetFlavourMC);

            flavourMap flavours;

			for (reco::JetFlavourMatchingCollection::const_iterator iFlavor = jetFlavourMC->begin(); iFlavor != jetFlavourMC->end(); iFlavor++) {
                flavours.insert(flavourMap::value_type(*((iFlavor->first).get()), abs(iFlavor->second.getFlavour())));
			}

			for (GenJetCollection::const_iterator jet_iter = GenJets->begin(); jet_iter!= GenJets->end(); ++jet_iter) {
				reco::GenJet myJet = reco::GenJet(*jet_iter);      
				if (myJet.pt() > 10) { 

					TCGenJet* jetCon = new ((*genJets)[genCount]) TCGenJet;
					jetCon->SetP4(myJet.px(), myJet.py(), myJet.pz(), myJet.energy());
					jetCon->SetHadEnergy(myJet.hadEnergy());
					jetCon->SetEmEnergy(myJet.emEnergy());
					jetCon->SetInvEnergy(myJet.invisibleEnergy());
					jetCon->SetAuxEnergy(myJet.auxiliaryEnergy());
					jetCon->SetNumConstit(myJet.getGenConstituents().size());

                    reco::Jet refJet(myJet.p4(),myJet.vertex());
                    if (flavours.find(refJet) != flavours.end()) jetCon->SetJetFlavor(flavours[refJet]);

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
			if (hlNames[i].compare(0, triggerPaths_[j].length(),triggerPaths_[j]) == 0) {
				triggerStatus |= 0x01 << j;
				pair<int, int> preScales;
				preScales = hltConfig_.prescaleValues(iEvent, iSetup, hlNames[i]); 
				hltPrescale[j] = preScales.first*preScales.second;
				//if (triggerPaths_[j] == "HLT_DoubleMu3_v") cout <<preScales.first<<"\t"<<preScales.second<<endl;
			}
		}
	} 

	++nEvents;

	if (eleCount > 0 || muCount > 0) eventTree -> Fill(); // possibly specify a cut in configuration

	primaryVtx->Clear("C");
	recoJets->Clear("C");
	recoMuons->Clear("C");
	recoElectrons->Clear("C");
	recoTaus->Clear("C");
	recoPhotons->Clear("C");
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
	recoElectrons  = new TClonesArray("TCElectron");
	recoMuons      = new TClonesArray("TCMuon");
	recoTaus       = new TClonesArray("TCTau");
	recoPhotons    = new TClonesArray("TCPhoton");
	genJets        = new TClonesArray("TCGenJet");
    genParticles   = new TClonesArray("TCGenParticle");
	beamSpot       = new TVector3();
	recoMET        = 0;

	eventTree->Branch("recoJets",&recoJets, 6400, 0);
	eventTree->Branch("recoElectrons",&recoElectrons, 6400, 0);
	eventTree->Branch("recoMuons",&recoMuons, 6400, 0);
	eventTree->Branch("recoTaus",&recoTaus, 6400, 0);
	eventTree->Branch("recoPhotons",&recoPhotons, 6400, 0);
	eventTree->Branch("recoMET", &recoMET, 6400, 0);
	eventTree->Branch("genJets",&genJets, 6400, 0);
	eventTree->Branch("genParticles",&genParticles, 6400, 0);

	eventTree->Branch("primaryVtx",&primaryVtx, 6400, 0);
	eventTree->Branch("beamSpot", &beamSpot, 6400, 0);
	eventTree->Branch("nPUVertices", &nPUVertices, "nPUVertices/I");

	eventTree->Branch("isRealData",&isRealData, "isRealData/O");
	eventTree->Branch("runNumber",&runNumber, "runNumber/I");
	eventTree->Branch("eventNumber",&eventNumber, "eventNumber/I");
	eventTree->Branch("lumiSection",&lumiSection, "lumiSection/I");
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
	eventTree->Branch("triggerStatus",&triggerStatus, "triggerStatus/i");
	eventTree->Branch("hltPrescale",hltPrescale, "hltPrescale[64]/i");

	runTree->Branch("deliveredLumi",&deliveredLumi, "deliveredLumi/F");
	runTree->Branch("recordedLumi",&recordedLumi, "recordedLumi/F");
	runTree->Branch("runNumber",&runNumber, "runNumber/i");

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
	edm::Handle<LumiSummary> lumiSummary;
	iLumi.getByLabel("lumiProducer", lumiSummary);

	deliveredLumi  += lumiSummary->avgInsDelLumi()*93.244;
	recordedLumi   += deliveredLumi*lumiSummary->liveFrac();
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

void ntupleProducer::associateJetToVertex(reco::PFJet inJet, Handle<reco::VertexCollection> vtxCollection, TCJet *outJet)
{
	const reco::TrackRefVector &tracks = inJet.getTrackRefs();

	vector<float>  associatedTrackSumPt;
	vector<const reco::Track*> jetTracks;
	float sumTrackX, sumTrackY, sumTrackZ, sumTrackPt;
	int   nJetTracks, nVertexTracks, nAssociatedTracks;
	int   vCount = 0;

	nJetTracks = nVertexTracks = nAssociatedTracks = 0;
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
		outJet->SetVtxIndex(0);
		outJet->SetVtx(0., 0., 0.);      	
    } else {
        outJet->SetVtx(sumTrackX/nJetTracks, sumTrackY/nJetTracks, sumTrackZ/nJetTracks);       
		for (VertexCollection::const_iterator iVtx = vtxCollection->begin(); iVtx!= vtxCollection->end(); ++iVtx) {	      
			reco::Vertex myVtx = reco::Vertex(*iVtx); 
			if(!myVtx.isValid() || myVtx.isFake()) continue;
			associatedTrackSumPt.push_back(0);            

			for(Vertex::trackRef_iterator iTrackRef = myVtx.tracks_begin(); iTrackRef != myVtx.tracks_end(); ++iTrackRef){
				const edm::RefToBase<reco::Track> &myTrackRef = *iTrackRef; 

				if(myTrackRef.isAvailable()){
					const reco::Track &myVertexTrack = *myTrackRef.get();		

					for(vector<const reco::Track*>::const_iterator iTrack = jetTracks.begin(); iTrack != jetTracks.end(); ++iTrack){
						if (*iTrack == &myVertexTrack) {
							associatedTrackSumPt.at(vCount) += myVertexTrack.pt()/sumTrackPt; 
							++nAssociatedTracks;
						}
					}
				}
			}
			++vCount;  
		}

		float maxSumPtFraction = 0;
		float maxCountFraction = 0;
		int   vertexIndex = 0;
		vCount = 0;

		for (vector<float>::const_iterator iTrackSumPt = associatedTrackSumPt.begin(); iTrackSumPt != associatedTrackSumPt.end(); ++iTrackSumPt) {
			if (*iTrackSumPt > maxSumPtFraction) {
				maxSumPtFraction = *iTrackSumPt;   
				vertexIndex      = vCount + 1;
			}
			++vCount;
		}
		outJet->SetVtxSumPtFrac(maxSumPtFraction);
		outJet->SetVtxSumPt(sumTrackPt);
		outJet->SetVtxTrackFrac((float)nAssociatedTracks/(float)nJetTracks);
		outJet->SetVtxNTracks(nJetTracks);
		outJet->SetVtxIndex(vertexIndex);
	}
}


bool ntupleProducer::isFilteredOutScraping( const edm::Event& iEvent, const edm::EventSetup& iSetup, int numtrack, double thresh)
{

	bool accepted = false;
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

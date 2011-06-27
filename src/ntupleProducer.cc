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
	electronIDMap_    = iConfig.getParameter<edm::InputTag>("electronIDMap");
	rhoCorrTag_       = iConfig.getUntrackedParameter<edm::InputTag>("rhoCorrTag");
	hlTriggerResults_ = iConfig.getUntrackedParameter<string>("HLTriggerResults","TriggerResults");
	hltProcess_       = iConfig.getUntrackedParameter<string>("hltName");
	triggerPaths_      = iConfig.getUntrackedParameter<vector<string> >("triggers");
	rootfilename      = iConfig.getUntrackedParameter<string>("rootfilename");

	saveJets_         = iConfig.getUntrackedParameter<bool>("saveJets");
	saveElectrons_    = iConfig.getUntrackedParameter<bool>("saveElectrons");
	saveMuons_        = iConfig.getUntrackedParameter<bool>("saveMuons");
	saveTaus_         = iConfig.getUntrackedParameter<bool>("saveTaus");
	savePhotons_      = iConfig.getUntrackedParameter<bool>("savePhotons");
	saveMET_          = iConfig.getUntrackedParameter<bool>("saveMET");
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

	int vtxCount, jetCount, metCount, muCount, eleCount, genCount, eleFakeCount, muFakeCount, partonCount;
	vtxCount = jetCount = metCount = muCount = eleCount = genCount = eleFakeCount = muFakeCount = partonCount = 0;
	float primaryVertexZ = -999;

	//////////////////////////
	//Get vertex information//
	//////////////////////////

	Handle<reco::VertexCollection> primaryVtcs;
	iEvent.getByLabel(primaryVtxTag_, primaryVtcs);

	for(VertexCollection::const_iterator vtx_iter = primaryVtcs->begin(); vtx_iter!= primaryVtcs->end(); ++vtx_iter){
		reco::Vertex myVtx = reco::Vertex(*vtx_iter);
		if(!myVtx.isValid() || myVtx.isFake()) continue;
		TCPrimaryVtx* vtxCon = new ((*primaryVtx)[vtxCount]) TCPrimaryVtx;
		vtxCon->SetPosition(myVtx.x(), myVtx.y(), myVtx.z());
		vtxCon->SetNDof(myVtx.ndof());
		vtxCon->SetChi2(myVtx.chi2());
		vtxCon->SetNTrks(myVtx.tracksSize()); 
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

		edm::Handle<reco::JetTagCollection> bTagHandle1;
		iEvent.getByLabel("trackCountingHighEffBJetTags", bTagHandle1);
		const reco::JetTagCollection & bTags1 = *(bTagHandle1.product());
		reco::JetTagCollection::const_iterator jet_it_1;

		//edm::ESHandle<JetCorrectorParametersCollection> JetCorParColl;
		//iSetup.get<JetCorrectionsRecord>().get("AK5PF",JetCorParColl);
		//JetCorrectorParameters const & JetCorPar = (*JetCorParColl)["Uncertainty"];
		//JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty(JetCorPar);

		const JetCorrector* correctorL1  = JetCorrector::getJetCorrector("ak5PFL1Fastjet",iSetup);
		const JetCorrector* correctorL2  = JetCorrector::getJetCorrector("ak5PFL2Relative",iSetup);
		const JetCorrector* correctorL3  = JetCorrector::getJetCorrector("ak5PFL3Absolute",iSetup);
		//const JetCorrector* correctorRes = JetCorrector::getJetCorrector("ak5PFResidual", iSetup);

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

			if (myJet.pt() < 10.) continue;
			//cout<<scale1<<", "<<scale2<<", "<<scale3<<"\n"<<endl;

			TCJet* jetCon = new ((*recoJets)[jetCount]) TCJet;

			jetCon->SetP4(myJet.px(), myJet.py(), myJet.pz(), myJet.energy());
			jetCon->SetVtx(-999.0, -999.0, -999.0);
			jetCon->SetChHadFrac(myJet.chargedHadronEnergyFraction());
			jetCon->SetNeuHadFrac(myJet.neutralHadronEnergyFraction());
			jetCon->SetChEmFrac(myJet.chargedEmEnergyFraction());
			jetCon->SetNeuEmFrac(myJet.neutralEmEnergyFraction());
			jetCon->SetNumConstit(myJet.chargedMultiplicity() + myJet.neutralMultiplicity());
			jetCon->SetNumChPart(myJet.chargedMultiplicity());

			// Get b-tag information
			for (jet_it_1 = bTags1.begin(); jet_it_1 != bTags1.end(); jet_it_1++) {
				if (sqrt(pow(jet_it_1->first->eta() - corJet.eta(), 2) + pow(deltaPhi(jet_it_1->first->phi(),corJet.phi()), 2)) == 0.) {
					jetCon->SetBDiscrTrkCountHiEff(jet_it_1->second);
				}
			}

			//cout<< "\t" << corJet.pt() << " | " << corJet.eta() << " | " << corJet.phi() <<endl;

			//add more corrections

			jetCon->SetJetCorr(1, scale1);
			jetCon->SetJetCorr(2, scale2);
			jetCon->SetJetCorr(3, scale3);

			//if (isRealData) {
			//	float scaleRes = correctorRes->correction(corJet);
			//	jetCon->SetJetCorr(4, scaleRes);
			//}

			//jecUnc->setJetEta(corJet.eta());
			//jecUnc->setJetPt(corJet.pt());
			//jetCon->SetUncertaintyJES(jecUnc->getUncertainty(true)); 

			/////////////////////////
			//get associated tracks//
			/////////////////////////


			const reco::TrackRefVector &tracks = myJet.getTrackRefs();

			vector<TVector3> vtxPositionCollection;
			vector<float>  associatedTrackSumPt;
			vector<const reco::Track*> jetTrackAddresses;
			float sumTrackX, sumTrackY, sumTrackZ, sumTrackIP, sumTrackPt;
			int   nJetTracks, nVertexTracks, nAssociatedTracks, vertexIndex;
			int   vCount = 0;

			nJetTracks = nVertexTracks = nAssociatedTracks = 0;
			sumTrackX = sumTrackY = sumTrackZ  = sumTrackIP  = sumTrackPt = 0;


			if(fabs(myJet.eta()) < 2.5){

				for (TrackRefVector::const_iterator iTrack = tracks.begin(); iTrack != tracks.end(); ++iTrack) {
					const reco::Track &myJetTrack = **iTrack;

					sumTrackPt += myJetTrack.pt();
					sumTrackX  += myJetTrack.vx();
					sumTrackY  += myJetTrack.vy();            
					sumTrackZ  += myJetTrack.vz();
					sumTrackIP += myJetTrack.dxy(vertexBeamSpot.position());
					jetTrackAddresses.push_back(&myJetTrack);
					++nJetTracks;
				}

				if (nJetTracks > 0) {
					jetCon->SetVtx(sumTrackX/nJetTracks, sumTrackY/nJetTracks, sumTrackZ/nJetTracks);      	
				}
				if(jetTrackAddresses.size() > 0){

					for (VertexCollection::const_iterator vtx_iter = primaryVtcs->begin(); vtx_iter!= primaryVtcs->end(); ++vtx_iter) {	      
						reco::Vertex myVtx = reco::Vertex(*vtx_iter); 
						if(!myVtx.isValid() || myVtx.isFake()) continue;
						TVector3 *iVtxPosition = new TVector3(myVtx.x(), myVtx.y(), myVtx.z());
						vtxPositionCollection.push_back(*iVtxPosition);
						associatedTrackSumPt.push_back(0);            
						for(Vertex::trackRef_iterator iTrackRef = myVtx.tracks_begin(); iTrackRef != myVtx.tracks_end(); ++iTrackRef){
							const edm::RefToBase<reco::Track> &myTrackRef = *iTrackRef; 

							if(myTrackRef.isAvailable()){
								const reco::Track &myVertexTrack = *myTrackRef.get();		

								for(vector<const reco::Track*>::const_iterator iTrackAddress = jetTrackAddresses.begin(); iTrackAddress != jetTrackAddresses.end(); ++iTrackAddress){
									if (*iTrackAddress == &myVertexTrack) {
										associatedTrackSumPt.at(vCount) += myVertexTrack.pt()/sumTrackPt; 
										++nAssociatedTracks;
									}
								}
							}
						}
						++vCount;  
					}

					float maxSumPtFraction = 0;
					vCount = vertexIndex = 0;

					for (vector<float>::const_iterator iTrackSumPt = associatedTrackSumPt.begin(); iTrackSumPt != associatedTrackSumPt.end(); ++iTrackSumPt) {
						if (*iTrackSumPt > maxSumPtFraction) {
							maxSumPtFraction = *iTrackSumPt;   
							vertexIndex      = vCount + 1;
						}
						++vCount;
					}
					jetCon->SetVtxSumPtFrac(maxSumPtFraction);
					jetCon->SetVtxSumPt(sumTrackPt);
					jetCon->SetVtxTrackFrac((float)nAssociatedTracks/(float)nJetTracks);
					jetCon->SetVtxNTracks(nJetTracks);
					jetCon->SetVtxIndex(vertexIndex);
				}
			} else {
				jetCon->SetVtxSumPtFrac(-1);
				jetCon->SetVtxSumPt(-1);
				jetCon->SetVtxTrackFrac(-1);
				jetCon->SetVtxNTracks(-1);
				jetCon->SetVtxIndex(0);
			}
			++jetCount;
		}   
	}



	/////////////
	// Get MET //
	/////////////

	if (saveMET_){

		Handle<PFMETCollection> MET;
		iEvent.getByLabel(metTag_, MET);

		if (pfMET.isValid()) {
			TCMET* metCon = new ((*recoMET)[metCount]) TCMET;
			metCon->SetSumEt(iMET->sumEt());
			metCon->SetMet(iMET->et());
			metCon->SetPhi(iMET->phi());
			metCon->SetPhotonEtFraction(iMET->photonEtFraction());
			metCon->SetElectronEtFraction(iMET->electronEtFraction());
			metCon->SetMuonEtFraction(iMET->muonEtFraction());
			metCon->SetNeutralHadronEtFraction(iMET->neutralHadronEtFraction());
			metCon->SetChargedHadronEtFraction(iMET->chargedHadronEtFraction());
			metCon->SetHFHadronEtFraction(iMET->HFHadronEtFraction());
			metCon->SetHFEMEtFraction(iMET->HFEMEtFraction());
			++metCount;
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
			muCon->Setp4(mu->px(), mu->py(), mu->pz(), mu->p());
			muCon->SetVtx(mu->globalTrack()->vx(),mu->globalTrack()->vy(),mu->globalTrack()->vz());
			muCon->SetCharge(mu->charge());
			muCon->SetisGLB(mu->isGlobalMuon());
			muCon->SetisTRK(mu->isTrackerMuon());
			muCon->Setdxy(mu->globalTrack()->dxy(vertexBeamSpot.position()));
			muCon->SetnPXLHits(mu->globalTrack()->hitPattern().numberOfValidPixelHits());
			muCon->SetnTRKHits(mu->globalTrack()->hitPattern().numberOfValidTrackerHits());
			muCon->SetnValidMuHits(mu->globalTrack()->hitPattern().numberOfValidMuonHits());
			muCon->SetNormChi2(mu->globalTrack()->normalizedChi2());
			muCon->SetnMatchSeg(mu->numberOfMatches());
			muCon->SetCaloComp(mu->caloCompatibility());
			muCon->SetSegComp(muon::segmentCompatibility(*mu));
			muCon->SetEMIso(mu->isolationR03().emEt);
			muCon->SetHADIso(mu->isolationR03().hadEt);
			muCon->SetTRKIso(mu->isolationR03().sumPt);
			muCount++;
		}
	}


	///////////////////
	// Get electrons //
	///////////////////

	if (saveElectrons_) {

		Handle<edm::ValueMap<float> > eIDValueMap;
		iEvent.getByLabel( electronIDMap_ , eIDValueMap );
		const edm::ValueMap<float> & eIDmap = * eIDValueMap ;

		Handle<GsfElectronCollection> electrons;
		iEvent.getByLabel(electronTag_, electrons);

		for (unsigned int i = 0; i < electrons->size(); i++){
			edm::Ref<reco::GsfElectronCollection> electronRef(electrons,i);
			//cout<<eIDmap[electronRef]<<"\t";

			if (electronRef->pt() < 10) continue;

			TCElectron* eleCon = new ((*recoElectrons)[eleCount]) TCElectron;
			eleCon->Setp4(electronRef->px(),electronRef->py(),electronRef->pz(),electronRef->p());
			eleCon->SetVtx(electronRef->gsfTrack()->vx(),electronRef->gsfTrack()->vy(),electronRef->gsfTrack()->vz());
			eleCon->SetCharge(electronRef->charge());
			eleCon->Setdxy(electronRef->gsfTrack()->dxy(vertexBeamSpot.position()));
			eleCon->SetNormChi2(electronRef->gsfTrack()->normalizedChi2());

			eleCon->SetIDMap(eIDmap[electronRef]);
			eleCon->SetEMIso(electronRef->dr03EcalRecHitSumEt());
			eleCon->SetHADIso(electronRef->dr03HcalTowerSumEt());
			eleCon->SetTRKIso(electronRef->dr03TkSumPt());
			eleCon->SetHoverE(electronRef->hadronicOverEm());
			eleCon->SetdPhiSC(electronRef->deltaPhiSuperClusterTrackAtVtx());
			eleCon->SetdEtaSC(electronRef->deltaEtaSuperClusterTrackAtVtx());
			eleCon->SetSig_IEtaIEta(electronRef->sigmaIetaIeta());
			eleCon->SetConversionFlag(electronRef->convFlags());
			eleCon->SetConversionDist(electronRef->convDist());
			eleCon->SetConversionDcot(electronRef->convDcot());
			eleCon->SetConversionRad(electronRef->convRadius());

			eleCon->SetPFChargedHadronIso(electronRef->pfIsolationVariables().chargedHadronIso);
			eleCon->SetPFNeutralHadronIso(electronRef->pfIsolationVariables().neutralHadronIso);
			eleCon->SetPFPhotonIso(electronRef->pfIsolationVariables().photonIso);

			eleCount++;
		}
	}

	/////////////////
	// Get photons //
	/////////////////


	//////////////
	// Get taus //
	//////////////


	////////////////////////
	// Get gen-level info //
	////////////////////////


	if (!isRealData) {

		Handle<HepMCProduct > genEvtHandle;
		iEvent.getByLabel( "generator", genEvtHandle) ;
		const HepMC::GenEvent* Evt = genEvtHandle->GetEvent();

		Handle<GenEventInfoProduct> GenEventInfoHandle;
		iEvent.getByLabel("generator", GenEventInfoHandle);

		Handle<reco::GenJetCollection> GenJets;
		iEvent.getByLabel(genJetTag_, GenJets);

		ptHat = qScale = -1; crossSection = 0;

		if (GenEventInfoHandle.isValid()) {
			qScale       = GenEventInfoHandle->qScale();
			ptHat        = (GenEventInfoHandle->hasBinningValues() ? GenEventInfoHandle->binningValues()[0] : 0.0);
			evtWeight    = GenEventInfoHandle->weight();

			h1_ptHat->Fill(ptHat);
		}

		vector<HepMC::GenParticle*> genPartons;

		for (HepMC::GenEvent::particle_const_iterator iGenParticle = Evt->particles_begin(); iGenParticle != Evt->particles_end(); ++iGenParticle) {
			HepMC::GenParticle *myGenPart = *iGenParticle;
			if (myGenPart->status() == 23) {
				new ((*hardPartonP4)[partonCount]) TLorentzVector(myGenPart->momentum().px(), myGenPart->momentum().py(), myGenPart->momentum().pz(), myGenPart->momentum().e());
				partonPdgId[partonCount] = myGenPart->pdg_id();
				++partonCount;
			}
		}
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

	for(int i = 0; i < (int)triggerPaths_.size(); ++i) 
		for (int i=0; i < (int)hlNames.size(); ++i) {      
			if (!triggerDecision(hltR, i)) continue;	
			for (int j = 0; j < (int)triggerPaths_.size(); ++j){
				hltPrescale[i] = hltConfig_.prescaleValue(iEvent, iSetup, triggerPaths_[i]); //This should be done at the end of the run
				if (hlNames[i].compare(0, triggerPaths_[j].length(),triggerPaths_[j]) == 0) {
					triggerStatus |= 0x01 << j;
				}
			}
		} 

	if (true) sTree -> Fill(); // possibly specify a cut in configuration

	primaryVtx->Clear("C");
	recoJets->Clear("C");
	genJets->Clear("C");
	recoMET->Clear("C");
	recoElectrons->Clear("C");
	recoMuons->Clear("C");
	hardPartonP4->Clear("C");

}

// ------------ method called once each job just before starting event loop  ------------
void  ntupleProducer::beginJob()
{  
	ntupleFile               = new TFile(rootfilename.c_str(), "RECREATE");
	sTree                    = new TTree("ntupleTree", "Tree for Jets");

	recoJets                 = new TClonesArray("TCJet");
	genJets                  = new TClonesArray("TCGenJet");
	recoMET                  = new TClonesArray("TCMET");
	recoElectrons            = new TClonesArray("TCElectron");
	recoMuons                = new TClonesArray("TCMuon");
	primaryVtx               = new TClonesArray("TCPrimaryVtx");
	hardPartonP4             = new TClonesArray("TLorentzVector");

	sTree->Branch("recoJets",&recoJets, 6400, 0);
	sTree->Branch("recoElectrons",&recoElectrons, 6400, 0);
	sTree->Branch("recoMuons",&recoMuons, 6400, 0);
	sTree->Branch("recoMET",&recoMET, 6400, 0);
	sTree->Branch("primaryVtx",&primaryVtx, 6400, 0);

	sTree->Branch("eventNumber",&eventNumber, "eventNumber/I");
	sTree->Branch("runNumber",&runNumber, "runNumber/I");
	sTree->Branch("lumiSection",&lumiSection, "lumiSection/I");
	sTree->Branch("triggerStatus",&triggerStatus, "triggerStatus/i");
	sTree->Branch("hltPrescale",hltPrescale, "hltPrescale[32]/i");
	sTree->Branch("isRealData",&isRealData, "isRealData/i");
	sTree->Branch("bunchCross",&bunchCross, "bunchCross/i");
	sTree->Branch("lumiDeadCount",&lumiDeadCount, "lumiDeadCount/f");
	sTree->Branch("lumiLiveFrac",&lumiLiveFrac, "lumiLiveFrac/f");
	sTree->Branch("intDelLumi",&intDelLumi, "intDelLumi/f");
	sTree->Branch("ptHat",&ptHat, "ptHat/f");
	sTree->Branch("qScale", &qScale, "qScale/f");
	sTree->Branch("evtWeight", &evtWeight, "evtWeight/f");
	sTree->Branch("crossSection", &crossSection, "crossSection/f");
	sTree->Branch("rhoFactor",&rhoFactor, "rhoFactor/F");
}

void ntupleProducer::beginRun(const edm::Run& iRun, const edm::EventSetup& iEvent)
{
	bool changed = true; 
	hltConfig_.init(iRun, iEvent, hltProcess_, changed);
}

void ntupleProducer::endLuminosityBlock(const edm::LuminosityBlock& iLumi, const edm::EventSetup& iEvent)
{
	edm::Handle<LumiSummary> lumiSummary;
	iLumi.getByLabel("lumiProducer", lumiSummary);

	lumiDeadCount  = lumiSummary->deadcount();
	lumiLiveFrac   = lumiSummary->liveFrac();
	intDelLumi     = lumiSummary->avgInsDelLumi()*93.244;

	//cout<<iLumi.id().luminosityBlock()<<endl;
	//cout<<"\t Dead Count = "<<lumiSummary->deadcount()<<endl;
	//cout<<"\t Fraction of dead time = "<<1 - lumiSummary->liveFrac()<<endl;
	//cout<<"\t Integrated luminosity = "<<lumiSummary->avgInsDelLumi()*93.244<<endl;
	//cout<<"\t Dead time corrected luminosity = "<<lumiSummary->avgInsDelLumi()*lumiSummary->liveFrac()*93.244<<endl;
}

void ntupleProducer::endRun(const edm::Run& iRun, const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	Handle<GenRunInfoProduct> GenRunInfoHandle;
	iEvent.getByLabel("generator", GenRunInfoHandle);

	if (GenRunInfoHandle.isValid()) {
		crossSection = GenRunInfoHandle->crossSection();
	}

}
// ------------ method called once each job just after ending the event loop  ------------
void ntupleProducer::endJob() 
{
	ntupleFile->Write();
	ntupleFile->Close();
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

//define this as a plug-in
DEFINE_FWK_MODULE(ntupleProducer);

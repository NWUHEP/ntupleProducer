// $Id: ntupleProducer.cc,v 1.9 2011/06/30 12:44:45 andrey Exp $

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
	saveGenJets_      = iConfig.getUntrackedParameter<bool>("saveGenJets");

	electronIDMap95_            = iConfig.getParameter<edm::InputTag>("electronIDMap95");
	electronIDMap90_            = iConfig.getParameter<edm::InputTag>("electronIDMap90");
	electronIDMap85_            = iConfig.getParameter<edm::InputTag>("electronIDMap85");
	electronIDMap80_            = iConfig.getParameter<edm::InputTag>("electronIDMap80");
	electronIDMap70_            = iConfig.getParameter<edm::InputTag>("electronIDMap70");
	electronIDMap60_            = iConfig.getParameter<edm::InputTag>("electronIDMap60");

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
		vtxCon->SetIsFake(myVtx.isFake());
		vtxCon->SetNtracks(myVtx.nTracks()); //0 is the minWeight, default is 0.5
		//if (myVtx.nTracks(0)!=myVtx.tracksSize()) LogWarning(" Vertex nTracks: ")<<"is something wrong here?";
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
		PFMETCollection::const_iterator pfMET = MET->begin();

		if (MET.isValid()) {
			TCMET* metCon = new ((*recoMET)[metCount]) TCMET;
			metCon->SetSumEt(pfMET->sumEt());
			metCon->SetMet(pfMET->et());
			metCon->SetPhi(pfMET->phi());
			metCon->SetPhotonEtFraction(pfMET->photonEtFraction());
			metCon->SetElectronEtFraction(pfMET->electronEtFraction());
			metCon->SetMuonEtFraction(pfMET->muonEtFraction());
			metCon->SetNeutralHadronEtFraction(pfMET->neutralHadronEtFraction());
			metCon->SetChargedHadronEtFraction(pfMET->chargedHadronEtFraction());
			metCon->SetHFHadronEtFraction(pfMET->HFHadronEtFraction());
			metCon->SetHFEMEtFraction(pfMET->HFEMEtFraction());
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

	  //Handle<edm::ValueMap<float> > eIDValueMap;
	  //iEvent.getByLabel( electronIDMap_ , eIDValueMap );
	  //const edm::ValueMap<float> & eIDmap = * eIDValueMap ;
	  
	  Handle<GsfElectronCollection> electrons;
	  iEvent.getByLabel(electronTag_, electrons);
	  
	  edm::Handle<edm::ValueMap<float> > eIDValueMap95;
	  iEvent.getByLabel( electronIDMap95_ , eIDValueMap95 );
	  const edm::ValueMap<float> & eIDmap95 = * eIDValueMap95 ;

	  edm::Handle<edm::ValueMap<float> > eIDValueMap90;
	  iEvent.getByLabel( electronIDMap90_ , eIDValueMap90 );
	  const edm::ValueMap<float> & eIDmap90 = * eIDValueMap90 ;
	  
	  edm::Handle<edm::ValueMap<float> > eIDValueMap85;
	  iEvent.getByLabel( electronIDMap85_ , eIDValueMap85 );
	  const edm::ValueMap<float> & eIDmap85 = * eIDValueMap85 ;
	  
	  edm::Handle<edm::ValueMap<float> > eIDValueMap80;
	  iEvent.getByLabel( electronIDMap80_ , eIDValueMap80 );
	  const edm::ValueMap<float> & eIDmap80 = * eIDValueMap80 ;
	  
	  edm::Handle<edm::ValueMap<float> > eIDValueMap70;
	  iEvent.getByLabel( electronIDMap70_ , eIDValueMap70 );
	  const edm::ValueMap<float> & eIDmap70 = * eIDValueMap70 ;
	  
	  edm::Handle<edm::ValueMap<float> > eIDValueMap60;
	  iEvent.getByLabel( electronIDMap60_ , eIDValueMap60 );
	  const edm::ValueMap<float> & eIDmap60 = * eIDValueMap60 ;
	  
	  
	  for (unsigned int i = 0; i < electrons->size(); i++){
	    edm::Ref<reco::GsfElectronCollection> electronRef(electrons,i);
	    //cout<<eIDmap[electronRef]<<"\t";
	    
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
	    
	    float sumPt3 = 0;
	    float gamma3 = 0;
	    float neutral3 = 0;
	    float sumPt4 = 0;
	    float gamma4 = 0;
	    float neutral4 = 0;
	    float sumPt5 = 0;
	    float gamma5 = 0;
	    float neutral5 = 0;
	    eleCon->SetPfSumPt(0.5, sumPt5);
	    eleCon->SetPfSumPt(0.4, sumPt4);
	    eleCon->SetPfSumPt(0.3, sumPt3);
	    eleCon->SetPfEGamma(0.5, gamma5);
	    eleCon->SetPfEGamma(0.4, gamma4);
	    eleCon->SetPfEGamma(0.3, gamma3);
	    eleCon->SetPfENeutral(0.5, neutral5);
	    eleCon->SetPfENeutral(0.4, neutral4);
	    eleCon->SetPfENeutral(0.3, neutral3);
	    
	    //eleCon->SetPfChargedHadronIso(electronRef->pfIsolationVariables().chargedHadronIso);
	    //eleCon->SetPfNeutralHadronIso(electronRef->pfIsolationVariables().neutralHadronIso);
	    //eleCon->SetPfPhotonIso(electronRef->pfIsolationVariables().photonIso);


	    eleCount++;
	  }
	}
	
	/////////////////
	// Get photons //
	/////////////////

	if (savePhotons_) {
		edm::Handle<reco::PhotonCollection> photons;
		iEvent.getByLabel(photonTag_, photons);

		int photonCount = 0;
		for (reco::PhotonCollection::const_iterator iPhoton = photons->begin(); iPhoton != photons->end() ; ++iPhoton)
		{
			TCPhoton* myPhoton = new ((*recoPhotons)[photonCount]) TCPhoton;
			myPhoton->SetP4(iPhoton->px(), iPhoton->py(), iPhoton->pz(), iPhoton->p());
			myPhoton->SetEMIso(iPhoton->ecalRecHitSumEtConeDR04());
			myPhoton->SetHADIso(iPhoton->hcalTowerSumEtConeDR04());
			myPhoton->SetTRKIso(iPhoton->trkSumPtHollowConeDR04());
			myPhoton->SetHadOverEm(iPhoton->hadronicOverEm());
			myPhoton->SetSigmaIEtaIEta(iPhoton->sigmaIetaIeta());
			//myPhoton->SetSigmaIphiIphi();
			//myPhoton->SetE2OverE9();
			myPhoton->SetEtaSupercluster(iPhoton->superCluster()->eta());
			myPhoton->SetTrackVeto(iPhoton->hasPixelSeed());

			photonCount++;
		}
	}

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

		evtWeight = ptHat = qScale = -1;

		if (GenEventInfoHandle.isValid()) {
			qScale       = GenEventInfoHandle->qScale();
			ptHat        = (GenEventInfoHandle->hasBinningValues() ? GenEventInfoHandle->binningValues()[0] : 0.0);
			//evtWeight    = GenEventInfoHandle->weight();
		}

		//////////////////////
		// Get genParticles //
		//////////////////////

		vector<HepMC::GenParticle*> genPartons;

		for (HepMC::GenEvent::particle_const_iterator iGenParticle = Evt->particles_begin(); iGenParticle != Evt->particles_end(); ++iGenParticle) {
			HepMC::GenParticle *myGenPart = *iGenParticle;
			if (myGenPart->status() == 23) {
				new ((*hardPartonP4)[partonCount]) TLorentzVector(myGenPart->momentum().px(), myGenPart->momentum().py(), myGenPart->momentum().pz(), myGenPart->momentum().e());
				partonPdgId[partonCount] = myGenPart->pdg_id();
				++partonCount;
			}
		}

		/////////////////
		// Get genJets //
		/////////////////

		if (saveGenJets_) {

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

					++genCount;	
				}
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

	for (int i=0; i < (int)hlNames.size(); ++i) {      
		if (!triggerDecision(hltR, i)) continue;	
		for (int j = 0; j < (int)triggerPaths_.size(); ++j){
			if (hlNames[i].compare(0, triggerPaths_[j].length(),triggerPaths_[j]) == 0) {
				triggerStatus |= 0x01 << j;
				pair<int, int> preScales;
				preScales = hltConfig_.prescaleValues(iEvent, iSetup, hlNames[i]); 
				hltPrescale[j] = preScales.first*preScales.second;
				cout<<hlNames[i]<<"\t"<<hltPrescale[j]<<endl;
			}
		}
	} 

	if (triggerStatus != 0x0) eventTree -> Fill(); // possibly specify a cut in configuration

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
	eventTree                = new TTree("eventTree", "Tree ");
	runTree                  = new TTree("runTree", "Tree for Jets");

	primaryVtx               = new TClonesArray("TCPrimaryVtx");
	recoJets                 = new TClonesArray("TCJet");
	recoMET                  = new TClonesArray("TCMET");
	recoElectrons            = new TClonesArray("TCElectron");
	recoMuons                = new TClonesArray("TCMuon");
	recoTaus                 = new TClonesArray("TCTau");
	recoPhotons              = new TClonesArray("TCPhoton");
	genJets                  = new TClonesArray("TCGenJet");
	hardPartonP4             = new TClonesArray("TLorentzVector");

	eventTree->Branch("primaryVtx",&primaryVtx, 6400, 0);
	eventTree->Branch("recoJets",&recoJets, 6400, 0);
	eventTree->Branch("recoElectrons",&recoElectrons, 6400, 0);
	eventTree->Branch("recoMuons",&recoMuons, 6400, 0);
	eventTree->Branch("recoTaus",&recoTaus, 6400, 0);
	eventTree->Branch("recoPhotons",&recoPhotons, 6400, 0);
	eventTree->Branch("recoMET",&recoMET, 6400, 0);
	eventTree->Branch("genJets",&genJets, 6400, 0);

	eventTree->Branch("eventNumber",&eventNumber, "eventNumber/I");
	eventTree->Branch("runNumber",&runNumber, "runNumber/I");
	eventTree->Branch("lumiSection",&lumiSection, "lumiSection/I");
	eventTree->Branch("triggerStatus",&triggerStatus, "triggerStatus/i");
	eventTree->Branch("isRealData",&isRealData, "isRealData/i");
	eventTree->Branch("bunchCross",&bunchCross, "bunchCross/i");
	eventTree->Branch("ptHat",&ptHat, "ptHat/f");
	eventTree->Branch("qScale", &qScale, "qScale/f");
	eventTree->Branch("evtWeight", &evtWeight, "evtWeight/f");
	eventTree->Branch("rhoFactor",&rhoFactor, "rhoFactor/F");
	eventTree->Branch("hltPrescale",hltPrescale, "hltPrescale[64]/i");

	runTree->Branch("deliveredLumi",&deliveredLumi, "deliveredLumi/f");
	runTree->Branch("recordedLumi",&recordedLumi, "recordedLumi/f");
	runTree->Branch("lumiDeadTime",&lumiDeadTime, "lumiDeadTime/f");
	runTree->Branch("runNumber",&runNumber, "runNumber/f");
	
	// Initialize HLT prescales //
	
	for (int i = 0; i < (int)triggerPaths_.size(); ++i) hltPrescale[i] = 1;
}

void ntupleProducer::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup)
{
	bool changed = true; 
	hltConfig_.init(iRun, iSetup, hltProcess_, changed);
	deliveredLumi = 0;
	recordedLumi  = 0;
	lumiDeadTime  = 0;
}

void ntupleProducer::endLuminosityBlock(const edm::LuminosityBlock& iLumi, const edm::EventSetup& iSetup)
{
	edm::Handle<LumiSummary> lumiSummary;
	iLumi.getByLabel("lumiProducer", lumiSummary);

	lumiDeadTime   += lumiSummary->deadcount()*93.244;
	deliveredLumi  += lumiSummary->avgInsDelLumi()*93.244;
	recordedLumi   += lumiSummary->avgInsDelLumi()*lumiSummary->liveFrac()*93.244;

	//cout<<iLumi.id().luminosityBlock()<<endl;
	//cout<<"\t Dead Count = "<<lumiSummary->deadcount()<<endl;
	//cout<<"\t Fraction of dead time = "<<1 - lumiSummary->liveFrac()<<endl;
	//cout<<"\t Integrated luminosity = "<<lumiSummary->avgInsDelLumi()*93.244<<endl;
	//cout<<"\t Dead time corrected luminosity = "<<lumiSummary->avgInsDelLumi()*lumiSummary->liveFrac()*93.244<<endl;
}

void ntupleProducer::endRun(const edm::Run& iRun, const edm::EventSetup& iSetup)
{
	cout<<"\t Integrated luminosity = "<<deliveredLumi<<endl;
	runTree->Fill();
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

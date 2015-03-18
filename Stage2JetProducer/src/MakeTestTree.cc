#include "Stage2Jets/Stage2JetProducer/interface/MakeTestTree.h"

#include "L1Trigger/L1TCalorimeter/plugins/TriggerTowerGeometry_new.cc"

#include <TMath.h>

//
// constructors and destructor
//
MakeTestTree::MakeTestTree(const edm::ParameterSet& iConfig)
{

  //  towersTag_ = iConfig.getParameter<edm::InputTag>("towerToken");
  m_towerToken = consumes<l1t::CaloTowerBxCollection>(towersTag_);

   //now do what ever initialization is needed
  edm::Service<TFileService> tfs;
  tree = tfs->make<TTree>("tree","tree");
  tree->Branch("run", &run);
  tree->Branch("event", &event);

  tree->Branch("towerPt",&towerPt);
  towerPt = new std::vector<int>();
  tree->Branch("towerEta",&towerEta);
  towerEta = new std::vector<int>();
  tree->Branch("towerPhi",&towerPhi);
  towerPhi = new std::vector<int>();
  tree->Branch("towerEtEm",&towerEtEm);
  towerEtEm = new std::vector<int>();
  tree->Branch("towerEtHad",&towerEtHad);
  towerEtHad = new std::vector<int>();

  tree->Branch("NjetUncalib",&NjetUncalib);
  jetUncalibPt = new std::vector<double>();
  tree->Branch("jetUncalibPt", &jetUncalibPt);
  jetUncalibPhi = new std::vector<double>();
  tree->Branch("jetUncalibPhi", &jetUncalibPhi);
  jetUncalibEta = new std::vector<double>();
  tree->Branch("jetUncalibEta", &jetUncalibEta);
  jetUncalibPhi_ph = new std::vector<double>();
  tree->Branch("jetUncalibPhi_ph", &jetUncalibPhi_ph);
  jetUncalibEta_ph = new std::vector<double>();
  tree->Branch("jetUncalibEta_ph", &jetUncalibEta_ph);


  tree->Branch("Njet",&Njet);
  jetPt = new std::vector<double>();
  tree->Branch("jetPt", &jetPt);
  jetPhi = new std::vector<double>();
  tree->Branch("jetPhi", &jetPhi);
  jetEta = new std::vector<double>();
  tree->Branch("jetEta", &jetEta);
  jetPhi_ph = new std::vector<double>();
  tree->Branch("jetPhi_ph", &jetPhi_ph);
  jetEta_ph = new std::vector<double>();
  tree->Branch("jetEta_ph", &jetEta_ph);

  tree->Branch("NjetNoPusUncalib",&NjetNoPusUncalib);
  jetNoPusUncalibPt = new std::vector<double>();
  tree->Branch("jetNoPusUncalibPt", &jetNoPusUncalibPt);
  jetNoPusUncalibPhi = new std::vector<double>();
  tree->Branch("jetNoPusUncalibPhi", &jetNoPusUncalibPhi);
  jetNoPusUncalibEta = new std::vector<double>();
  tree->Branch("jetNoPusUncalibEta", &jetNoPusUncalibEta);
  jetNoPusUncalibPhi_ph = new std::vector<double>();
  tree->Branch("jetNoPusUncalibPhi_ph", &jetNoPusUncalibPhi_ph);
  jetNoPusUncalibEta_ph = new std::vector<double>();
  tree->Branch("jetNoPusUncalibEta_ph", &jetNoPusUncalibEta_ph);

  tree->Branch("NjetNoPus",&NjetNoPus);
  jetNoPusPt = new std::vector<double>();
  tree->Branch("jetNoPusPt", &jetNoPusPt);
  jetNoPusPhi = new std::vector<double>();
  tree->Branch("jetNoPusPhi", &jetNoPusPhi);
  jetNoPusEta = new std::vector<double>();
  tree->Branch("jetNoPusEta", &jetNoPusEta);
  jetNoPusPhi_ph = new std::vector<double>();
  tree->Branch("jetNoPusPhi_ph", &jetNoPusPhi_ph);
  jetNoPusEta_ph = new std::vector<double>();
  tree->Branch("jetNoPusEta_ph", &jetNoPusEta_ph);

  tree->Branch("NjetGlobalUncalib",&NjetGlobalUncalib);
  jetGlobalUncalibPt = new std::vector<double>();
  tree->Branch("jetGlobalUncalibPt", &jetGlobalUncalibPt);
  jetGlobalUncalibPhi = new std::vector<double>();
  tree->Branch("jetGlobalUncalibPhi", &jetGlobalUncalibPhi);
  jetGlobalUncalibEta = new std::vector<double>();
  tree->Branch("jetGlobalUncalibEta", &jetGlobalUncalibEta);
  jetGlobalUncalibPhi_ph = new std::vector<double>();
  tree->Branch("jetGlobalUncalibPhi_ph", &jetGlobalUncalibPhi_ph);
  jetGlobalUncalibEta_ph = new std::vector<double>();
  tree->Branch("jetGlobalUncalibEta_ph", &jetGlobalUncalibEta_ph);

  tree->Branch("NjetGlobal",&NjetGlobal);
  jetGlobalPt = new std::vector<double>();
  tree->Branch("jetGlobalPt", &jetGlobalPt);
  jetGlobalPhi = new std::vector<double>();
  tree->Branch("jetGlobalPhi", &jetGlobalPhi);
  jetGlobalEta = new std::vector<double>();
  tree->Branch("jetGlobalEta", &jetGlobalEta);
  jetGlobalPhi_ph= new std::vector<double>();
  tree->Branch("jetGlobalPhi_ph", &jetGlobalPhi_ph);
  jetGlobalEta_ph = new std::vector<double>();
  tree->Branch("jetGlobalEta_ph", &jetGlobalEta_ph);

  tree->Branch("mht", &mht);
  tree->Branch("mhtPhi", &mhtPhi);
  tree->Branch("ht", &ht);

  tree->Branch("met", &met);
  tree->Branch("metPhi", &metPhi);
  tree->Branch("et", &et);

  tree->Branch("mhtNoPus", &mhtNoPus);
  tree->Branch("mhtPhiNoPus", &mhtPhiNoPus);
  tree->Branch("htNoPus", &htNoPus);

  tree->Branch("mhtGlobal", &mhtGlobal);
  tree->Branch("mhtPhiGlobal", &mhtPhiGlobal);
  tree->Branch("htGlobal", &htGlobal);

}


MakeTestTree::~MakeTestTree()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  delete towerPt;
  delete towerEta;
  delete towerPhi;
  delete towerEtEm;
  delete towerEtHad;

  delete jetPt;
  delete jetPhi;
  delete jetEta;
  delete jetUncalibPt;
  delete jetUncalibPhi;
  delete jetUncalibEta;
  delete jetNoPusPt;
  delete jetNoPusPhi;
  delete jetNoPusEta;
  delete jetNoPusUncalibPt;
  delete jetNoPusUncalibPhi;
  delete jetNoPusUncalibEta;
  delete jetGlobalPt;
  delete jetGlobalPhi;
  delete jetGlobalEta;
  delete jetGlobalUncalibPt;
  delete jetGlobalUncalibPhi;
  delete jetGlobalUncalibEta;
}


//
// member functions
//

// ------------ method called for each event  ------------
void
MakeTestTree::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
  towerPt->clear();
  towerEta->clear();
  towerPhi->clear();
  towerEtEm->clear();
  towerEtHad->clear();

  jetPt->clear();
  jetPhi->clear();
  jetEta->clear();
  jetPhi_ph->clear();
  jetEta_ph->clear();

  jetUncalibPt->clear();
  jetUncalibPhi->clear();
  jetUncalibEta->clear();
  jetUncalibPhi_ph->clear();
  jetUncalibEta_ph->clear();

  jetNoPusPt->clear();
  jetNoPusPhi->clear();
  jetNoPusEta->clear();
  jetNoPusPhi_ph->clear();
  jetNoPusEta_ph->clear();

  jetNoPusUncalibPt->clear();
  jetNoPusUncalibPhi->clear();
  jetNoPusUncalibEta->clear();
  jetNoPusUncalibPhi_ph->clear();
  jetNoPusUncalibEta_ph->clear();

  jetGlobalPt->clear();
  jetGlobalPhi->clear();
  jetGlobalEta->clear();
  jetGlobalPhi_ph->clear();
  jetGlobalEta_ph->clear();

  jetGlobalUncalibPt->clear();
  jetGlobalUncalibPhi->clear();
  jetGlobalUncalibEta->clear();
  jetGlobalUncalibPhi_ph->clear();
  jetGlobalUncalibEta_ph->clear();

  run = 0;
  event = 0;
  mht=-10;
  mhtPhi=-10;
  ht=-10;

   //Get the trigger towers and put them in an array
  // edm::Handle<l1slhc::L1CaloTowerCollection> triggerTowers;
  // iEvent.getByLabel("L1CaloTowerProducer", triggerTowers);


  edm::Handle< BXVector<l1t::CaloTower> > triggerTowers;
  iEvent.getByLabel("caloStage2Digis", "MP",triggerTowers);

  edm::Handle<std::vector<L1JetParticle> > jetHandle;
  iEvent.getByLabel("Stage2JetProducer","l1Stage2JetsDonutPUS", jetHandle);

  edm::Handle<std::vector<L1JetParticle> > jetUncalibHandle;
  iEvent.getByLabel("Stage2JetProducer","l1Stage2JetsDonutPUSUncalib", jetUncalibHandle);

  edm::Handle<std::vector<L1JetParticle> > jetGlobalHandle;
  iEvent.getByLabel("Stage2JetProducer","l1Stage2JetsGlobalPUS", jetGlobalHandle);

  edm::Handle<std::vector<L1JetParticle> > jetGlobalUncalibHandle;
  iEvent.getByLabel("Stage2JetProducer","l1Stage2JetsGlobalPUSUncalib", jetGlobalUncalibHandle);

  edm::Handle<std::vector<L1JetParticle> > jetNoPusHandle;
  iEvent.getByLabel("Stage2JetProducer","l1Stage2JetsNoPUS", jetNoPusHandle);

  edm::Handle<std::vector<L1JetParticle> > jetNoPusUncalibHandle;
  iEvent.getByLabel("Stage2JetProducer","l1Stage2JetsNoPUSUncalib", jetNoPusUncalibHandle);

  edm::Handle<std::vector<L1EtMissParticle> > mhtHandle;
  iEvent.getByLabel("Stage2JetProducer","l1Stage2DonutPUSMht", mhtHandle);

  edm::Handle<std::vector<L1EtMissParticle> > metHandle;
  iEvent.getByLabel("Stage2JetProducer","l1Stage2Met", metHandle);

  edm::Handle<std::vector<L1EtMissParticle> > mhtNoPusHandle;
  iEvent.getByLabel("Stage2JetProducer","l1Stage2NoPUSMht", mhtNoPusHandle);

  edm::Handle<std::vector<L1EtMissParticle> > mhtGlobalHandle;
  iEvent.getByLabel("Stage2JetProducer","l1Stage2GlobalPUSMht", mhtGlobalHandle);

  int nj=0;

  TriggerTowerGeometry ttg;

  double ET=0.;
  // int32 towPhi;
  // Trigger towers
  //  for(auto j=triggerTowers->begin(); j!=triggerTowers->end(); ++j) {
  for ( int ibx=triggerTowers->getFirstBX(); ibx<=triggerTowers->getLastBX(); ++ibx) {
      //  for(auto j=triggerTowers->begin(); j!=triggerTowers->end(); j++) {
     for ( auto j = triggerTowers->begin(ibx); j !=triggerTowers->end(ibx); ++j ) {

       if ( abs((*j).hwEta()) > 28 ) { continue; } 

       ET = ((*j).hwPt()); //(*j).E()+(*j).H();

       //  towPhi=(*j).hwPhi();

       // if(towPhi <= 0) { towPhi = towPhi+72; } 
       // else if (towPhi > 72) { towPhi = towPhi-72; } 
       // else { towPhi = towPhi; }


       // if (ET==0) continue;
       towerPt->push_back(ET);
       towerEta->push_back((*j).hwEta());
       towerPhi->push_back((*j).hwPhi()); //towPhi);
              

       // towerEtEm->push_back((*j).hwEtEm());
       // towerEtHad->push_back((*j).hwEtHad());
       
       nj++;
     }
  }

  int nMax=4;

  //  double physical_phi;
  std::cout << "" << std::endl;                                                                                    
  std::cout << "Calib chunky  L1Jets:" << std::endl;  
  
  nj=0;
  for (unsigned i = 0; i < jetHandle->size(); ++i) {
    
    L1JetParticle jet = jetHandle->at(i); 
    if (jet.pt()<10.) continue;

    if (nj>=nMax) break;
  
    jetPt->push_back(jet.pt());
    jetPhi->push_back(ttg.iPhi(jet.phi()));
    jetEta->push_back(ttg.iEta(jet.eta()));

    jetPhi_ph->push_back(jet.phi());
    jetEta_ph->push_back(jet.eta());

  
    std::cout << "Jet #" << nj << " has pT= " << jet.pt()  << 
      " , hw_eta = " << ttg.iEta(jet.eta()) << " , and hw_phi= " << ttg.iPhi(jet.phi()) << std::endl;
  
    nj++;
  }
  Njet=nj;

  std::cout << "" << std::endl;
  std::cout << "Uncalib chunky L1Jets:" << std::endl;
  nj=0;
  for (unsigned i = 0; i < jetUncalibHandle->size(); ++i) {

    L1JetParticle jetUncalib = jetUncalibHandle->at(i);
    if (jetUncalib.pt()<10.) continue;

    if (nj>=nMax) break;

    jetUncalibPt->push_back(jetUncalib.pt());
    jetUncalibPhi->push_back(ttg.iPhi(jetUncalib.phi()) );
    jetUncalibEta->push_back(ttg.iEta(jetUncalib.eta()));

    jetUncalibPhi_ph->push_back(jetUncalib.phi() );
    jetUncalibEta_ph->push_back(jetUncalib.eta());
  
    std::cout << "Jet #" << nj << " has pT= " << jetUncalib.pt()  << 
      " , hw_eta = " << ttg.iEta(jetUncalib.eta()) << " , and hw_phi= " << ttg.iPhi(jetUncalib.phi()) << std::endl; 
  
    // std::cout << "Jet #" << nj << " has pT= " << jetUncalib.pt()  << 
    //   " , eta = " << jetUncalib.eta() << " , and phi= " << physical_phi << std::endl;

    nj++;
  }
  NjetUncalib=nj;


  //  std::cout << "" << std::endl;
  // std::cout << "NoPusCalib L1Jets:" << std::endl; 
  nj=0;
  for (unsigned i = 0; i < jetNoPusHandle->size(); ++i) {
    L1JetParticle jetNoPus = jetNoPusHandle->at(i);

    if (jetNoPus.pt()<10.) continue;

    if (nj>=nMax) break;

    jetNoPusPt->push_back(jetNoPus.pt());
    jetNoPusPhi->push_back(2*ttg.iPhi(jetNoPus.phi()));
    jetNoPusEta->push_back(ttg.iEta(jetNoPus.eta()));

    jetNoPusPhi_ph->push_back(jetNoPus.phi());
    jetNoPusEta_ph->push_back(jetNoPus.eta());
    /*
    std::cout << "Jet #" << nj << " has pT= " << jetNoPus.pt()  << 
      " , hw_eta = " << ttg.iEta(jetNoPus.eta()) << " , and hw_phi= " << ttg.iPhi(jetNoPus.phi()) << std::endl;
    */
    nj++;
  }
  NjetNoPus=nj;


  std::cout << "" << std::endl;
  std::cout << "NoPusUncalib L1Jets:" << std::endl;
  nj=0;
  for (unsigned i = 0; i < jetNoPusUncalibHandle->size(); ++i) {

    L1JetParticle jetNoPusUncalib = jetNoPusUncalibHandle->at(i);

    if (jetNoPusUncalib.pt()<10.) continue;
    if (nj>=nMax) break;

    jetNoPusUncalibPt->push_back(jetNoPusUncalib.pt());
    jetNoPusUncalibPhi->push_back(2*ttg.iPhi(jetNoPusUncalib.phi()));
    jetNoPusUncalibEta->push_back(ttg.iEta(jetNoPusUncalib.eta()));

    jetNoPusUncalibPhi_ph->push_back(jetNoPusUncalib.phi());
    jetNoPusUncalibEta_ph->push_back(jetNoPusUncalib.eta());
 
    std::cout << "Jet #" << nj << " has pT= " << jetNoPusUncalib.pt()  << 
      " , hw_eta = " << ttg.iEta(jetNoPusUncalib.eta()) << " , and hw_phi= " << ttg.iPhi(jetNoPusUncalib.phi()) << std::endl;
 
    nj++;
  }
  NjetNoPusUncalib=nj;


  nj=0;
  for (unsigned i = 0; i < jetGlobalHandle->size(); ++i) {

    L1JetParticle jetGlobal = jetGlobalHandle->at(i);
    // if (jetGlobal.pt()<10.) continue;

    if (nj>=nMax) break;

    jetGlobalPt->push_back(jetGlobal.pt());
    jetGlobalPhi->push_back(2*ttg.iPhi(jetGlobal.phi()));
    jetGlobalEta->push_back(ttg.iEta(jetGlobal.eta()));

    jetGlobalPhi_ph->push_back(jetGlobal.phi());
    jetGlobalEta_ph->push_back(jetGlobal.eta());

    nj++;
  }
  NjetGlobal=nj;

  nj=0;
  for (unsigned i = 0; i < jetGlobalUncalibHandle->size(); ++i) {
    L1JetParticle jetGlobalUncalib = jetGlobalUncalibHandle->at(i);

    if (jetGlobalUncalib.pt()<10.)continue;
    if (nj>=nMax) break;

    jetGlobalUncalibPt->push_back(jetGlobalUncalib.pt());
    jetGlobalUncalibPhi->push_back(2*ttg.iPhi(jetGlobalUncalib.phi()));
    jetGlobalUncalibEta->push_back(ttg.iEta(jetGlobalUncalib.eta()));

    jetGlobalUncalibPhi_ph->push_back(jetGlobalUncalib.phi());
    jetGlobalUncalibEta_ph->push_back(jetGlobalUncalib.eta());

    nj++;
  }
  NjetGlobalUncalib=nj;

  // double d_mhtPhi, d_metPhi, d_mhtPhiNoPus, d_mhtPhiGlobal;

  mht=mhtHandle->at(0).pt();
  mhtPhi=ttg.iPhi(mhtHandle->at(0).phi());
  // if (d_mhtPhi<0) d_mhtPhi+=2.*TMath::Pi();
  // mhtPhi=d_mhtPhi;
  ht=mhtHandle->at(0).etTotal();

  met=metHandle->at(0).pt();
  metPhi=ttg.iPhi(metHandle->at(0).phi());
  // if (d_metPhi<0) d_metPhi+=2.*TMath::Pi();
  // metPhi=d_metPhi;

  // std::cout << "producer MET = " <<  met << " , and has PHI = " <<  metPhi << std::endl;
  // //	    " , and rounded physical PHI= " << metPhi << std::endl;
  // std::cout << "producer MHT = " <<  mht << " , and has PHI = " <<  mhtPhi << std::endl;

  et=metHandle->at(0).etTotal();

  mhtNoPus=mhtNoPusHandle->at(0).pt();
  mhtPhiNoPus=ttg.iPhi(mhtNoPusHandle->at(0).phi());
  // if (d_mhtPhiNoPus<0) d_mhtPhiNoPus+=2.*TMath::Pi();
  // mhtPhiNoPus=d_mhtPhiNoPus;

  htNoPus=mhtNoPusHandle->at(0).etTotal();
  
  mhtGlobal=mhtGlobalHandle->at(0).pt();
  mhtPhiGlobal=ttg.iPhi(mhtGlobalHandle->at(0).phi());
  // if (d_mhtPhiGlobal<0) d_mhtPhiGlobal+=2.*TMath::Pi();
  // mhtPhiGlobal=d_mhtPhiGlobal;

  htGlobal=mhtGlobalHandle->at(0).etTotal();

  tree->Fill();
}


// ------------ method called once each job just before starting event loop  ------------
  void 
MakeTestTree::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
  void 
MakeTestTree::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
  void 
MakeTestTree::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
  void 
MakeTestTree::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
  void 
MakeTestTree::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
  void 
MakeTestTree::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MakeTestTree::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MakeTestTree);

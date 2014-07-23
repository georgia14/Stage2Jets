#include "Stage2Jets/Stage2JetProducer/interface/MakeTestTree.h"

//
// constructors and destructor
//
MakeTestTree::MakeTestTree(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed
  edm::Service<TFileService> tfs;
  tree = tfs->make<TTree>("tree","tree");
  tree->Branch("run", &run);
  tree->Branch("event", &event);

  jetUncalibPt = new std::vector<double>();
  tree->Branch("jetUncalibPt", &jetUncalibPt);
  jetUncalibPhi = new std::vector<double>();
  tree->Branch("jetUncalibPhi", &jetUncalibPhi);
  jetUncalibEta = new std::vector<double>();
  tree->Branch("jetUncalibEta", &jetUncalibEta);

  jetPt = new std::vector<double>();
  tree->Branch("jetPt", &jetPt);
  jetPhi = new std::vector<double>();
  tree->Branch("jetPhi", &jetPhi);
  jetEta = new std::vector<double>();
  tree->Branch("jetEta", &jetEta);

  jetNoPusUncalibPt = new std::vector<double>();
  tree->Branch("jetNoPusUncalibPt", &jetNoPusUncalibPt);
  jetNoPusUncalibPhi = new std::vector<double>();
  tree->Branch("jetNoPusUncalibPhi", &jetNoPusUncalibPhi);
  jetNoPusUncalibEta = new std::vector<double>();
  tree->Branch("jetNoPusUncalibEta", &jetNoPusUncalibEta);

  jetNoPusPt = new std::vector<double>();
  tree->Branch("jetNoPusPt", &jetNoPusPt);
  jetNoPusPhi = new std::vector<double>();
  tree->Branch("jetNoPusPhi", &jetNoPusPhi);
  jetNoPusEta = new std::vector<double>();
  tree->Branch("jetNoPusEta", &jetNoPusEta);


  tree->Branch("mht", &mht);
  tree->Branch("mhtPhi", &mhtPhi);
  tree->Branch("ht", &ht);

  tree->Branch("mhtNoPus", &mhtNoPus);
  tree->Branch("mhtPhiNoPus", &mhtPhiNoPus);
  tree->Branch("htNoPus", &htNoPus);


}


MakeTestTree::~MakeTestTree()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
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
}


//
// member functions
//

// ------------ method called for each event  ------------
void
MakeTestTree::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  jetPt->clear();
  jetPhi->clear();
  jetEta->clear();
  jetUncalibPt->clear();
  jetUncalibPhi->clear();
  jetUncalibEta->clear();
  jetNoPusPt->clear();
  jetNoPusPhi->clear();
  jetNoPusEta->clear();
  jetNoPusUncalibPt->clear();
  jetNoPusUncalibPhi->clear();
  jetNoPusUncalibEta->clear();
  run = 0;
  event = 0;
  mht=-10;
  mhtPhi=-10;
  ht=-10;

  edm::Handle<std::vector<L1JetParticle> > jetHandle;
  iEvent.getByLabel("Stage2JetProducer","l1Stage2JetsDonutPUS", jetHandle);

  edm::Handle<std::vector<L1JetParticle> > jetUncalibHandle;
  iEvent.getByLabel("Stage2JetProducer","l1Stage2JetsDonutPUSUncalib", jetUncalibHandle);

  edm::Handle<std::vector<L1JetParticle> > jetNoPusHandle;
  iEvent.getByLabel("Stage2JetProducer","l1Stage2JetsNoPUS", jetNoPusHandle);

  edm::Handle<std::vector<L1JetParticle> > jetNoPusUncalibHandle;
  iEvent.getByLabel("Stage2JetProducer","l1Stage2JetsNoPUSUncalib", jetNoPusUncalibHandle);

  edm::Handle<std::vector<L1EtMissParticle> > mhtHandle;
  iEvent.getByLabel("Stage2JetProducer","l1Stage2DonutPUSMht", mhtHandle);

  edm::Handle<std::vector<L1EtMissParticle> > mhtNoPusHandle;
  iEvent.getByLabel("Stage2JetProducer","l1Stage2NoPUSMht", mhtNoPusHandle);
  

  for (unsigned i = 0; i < jetHandle->size(); ++i) {
    L1JetParticle jet = jetHandle->at(i);
    jetPt->push_back(jet.pt());
    jetPhi->push_back(jet.phi());
    jetEta->push_back(jet.eta());
  }

  for (unsigned i = 0; i < jetUncalibHandle->size(); ++i) {
    L1JetParticle jetUncalib = jetUncalibHandle->at(i);
    jetUncalibPt->push_back(jetUncalib.pt());
    jetUncalibPhi->push_back(jetUncalib.phi());
    jetUncalibEta->push_back(jetUncalib.eta());
  }

  for (unsigned i = 0; i < jetNoPusHandle->size(); ++i) {
    L1JetParticle jetNoPus = jetNoPusHandle->at(i);
    jetNoPusPt->push_back(jetNoPus.pt());
    jetNoPusPhi->push_back(jetNoPus.phi());
    jetNoPusEta->push_back(jetNoPus.eta());
  }

  for (unsigned i = 0; i < jetNoPusUncalibHandle->size(); ++i) {
    L1JetParticle jetNoPusUncalib = jetNoPusUncalibHandle->at(i);
    jetNoPusUncalibPt->push_back(jetNoPusUncalib.pt());
    jetNoPusUncalibPhi->push_back(jetNoPusUncalib.phi());
    jetNoPusUncalibEta->push_back(jetNoPusUncalib.eta());
  }

  mht=mhtHandle->at(0).pt();
  mhtPhi=mhtHandle->at(0).phi();

  ht=mhtHandle->at(0).etTotal();

  mhtNoPus=mhtNoPusHandle->at(0).pt();
  mhtPhiNoPus=mhtNoPusHandle->at(0).phi();

  htNoPus=mhtNoPusHandle->at(0).etTotal();

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

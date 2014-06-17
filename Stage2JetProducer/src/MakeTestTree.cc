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
  jetPt = new std::vector<double>();
  tree->Branch("jetPt", &jetPt);
  jetPhi = new std::vector<double>();
  tree->Branch("jetPhi", &jetPhi);
  jetEta = new std::vector<double>();
  tree->Branch("jetEta", &jetEta);
  tree->Branch("mht", &mht);
  tree->Branch("mhtPhi", &mhtPhi);
  tree->Branch("ht", &ht);
  

}


MakeTestTree::~MakeTestTree()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
  delete jetPt;
  delete jetPhi;
  delete jetEta;
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
  run = 0;
  event = 0;
  mht=-10;
  mhtPhi=-10;
  ht=-10;

  edm::Handle<std::vector<L1JetParticle> > jetHandle;
  iEvent.getByLabel("Stage2JetProducer","l1Stage2Jets", jetHandle);

  edm::Handle<std::vector<L1EtMissParticle> > mhtHandle;
  iEvent.getByLabel("Stage2JetProducer","l1Stage2Mht", mhtHandle);

  edm::Handle<double> htHandle;
  iEvent.getByLabel("Stage2JetProducer","l1Stage2Ht", htHandle);

  for (unsigned i = 0; i < jetHandle->size(); ++i) {
    L1JetParticle jet = jetHandle->at(i);
    jetPt->push_back(jet.pt());
    jetPhi->push_back(jet.phi());
    jetEta->push_back(jet.eta());
  }

  mht=mhtHandle->at(0).pt();
  mhtPhi=mhtHandle->at(0).phi();

  ht=*htHandle;

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

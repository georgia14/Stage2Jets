#include "Stage2Jets/Stage2JetProducer/interface/Stage2JetProducer.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

namespace{
  bool sortbypt(const LeafCandidate &a, const LeafCandidate &b) { return a.pt() > b.pt(); }
}

Stage2JetProducer::Stage2JetProducer(const edm::ParameterSet& iConfig) {

  //Produces the jets and MET etc
  produces<std::vector<L1JetParticle> >("l1Stage2Jets");
  produces<std::vector<L1JetParticle> >("l1Stage2JetsUncalib");
  produces<std::vector<reco::LeafCandidate> >("l1Stage2Mht");
  produces<double>("l1Stage2Ht");

  //Get the configs
  towersTag_ = iConfig.getParameter<edm::InputTag>("towerInput");
  mhtThreshold_ = iConfig.getParameter<double>("mhtThreshold");
  htThreshold_ = iConfig.getParameter<double>("htThreshold");

}


Stage2JetProducer::~Stage2JetProducer() {
}

// ------------ method called to produce the data  ------------
void Stage2JetProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  //Get the trigger towers and put them in an array
  edm::Handle<l1slhc::L1CaloTowerCollection> triggerTowers;
  iEvent.getByLabel(towersTag_, triggerTowers);

  //Maybe make these configurable but not for now
  int jetsize=5;
  int vetowindowsize=4;
  int seedthresh1=0;
  int seedthresh2=0;

  //Make these incase you want them in the future
  double ET=0.;
  double met_x=0.;
  double met_y=0.;

  std::vector< std::vector<int> > ttArray(56, std::vector<int>(72, 0)); //this is just a container for the (E+H) per tower

  //std::auto_ptr<std::vector<reco::LeafCandidate> > uncalibL1Jets(new
  //    std::vector<reco::LeafCandidate>());
  std::vector<L1JetParticle> uncalibL1Jets;

  for(auto j=triggerTowers->begin(); j!=triggerTowers->end(); j++) {

    if ( abs((*j).iEta()) > 28 ) { continue; } //i.e. |eta| < 3 only

    ttArray[g.new_iEta((*j).iEta())][g.new_iPhi((*j).iPhi())] = ((*j).E() + (*j).H());
    ET += (*j).E()+(*j).H();
    met_x -= cos(g.phi((*j).iPhi())) * ((*j).E() + (*j).H());
    met_y -= sin(g.phi((*j).iPhi())) * ((*j).E() + (*j).H());
    //so now ttArray is on the scale ieta 0-56, iphi 0-71
  }
 


  std::vector<int> areas(jetsize+1,0); //to hold the ring areas (i.e. when we get up against the boundaries)
  std::vector<int> jetareas; //to hold the ring areas (i.e. when we get up against the boundaries)

  for ( int i = 0; i < (int)ttArray.size(); i++) {
    for ( int j = 0; j < (int)ttArray[i].size(); j++) {
      std::vector<int> jetTower;
      //std::cout << "new: (" << i << ", " << j << ", " << ttArray[i][j] << ")" << std::endl;
      int numtowersaboveme=0;
      int numtowersabovethresh=0;
      //int seedtower = ttArray[i][j];  
      std::vector<int> localsums(jetsize+1,0); //to hold the ring sums (+1 for centre)
      std::vector<std::pair<int,int>> outerstrips(4,std::make_pair(0,0)); //to hold the energies in the 4 surrounding outer strips (excluding corners)
      int jetarea = 1;
      //int pusarea = 0;
      for(int k=(i-jetsize); k<=(i+jetsize); k++) {
        for(int l=(j-jetsize); l<=(j+jetsize); l++) {
          if(k < 0 || k > 55) continue; //i.e. out of bounds of eta<3
          //std::cout << " k = " << k << ", l = " << l << ", i =" << i << ", j = " << j << std::endl;

          //make a co-ordinate transform at the phi boundary
          int newl;
          if(l < 0) { newl = l+72; } 
          else if (l > 71) { newl = l-72; } 
          else { newl = l; }
          if (l != j && k != i)
          {
            jetTower.push_back(ttArray[k][newl]);
          }
          if(ttArray[k][newl] > seedthresh2) { numtowersabovethresh++; }

          //to decide which ring to assign energy to
          for( int m=0; m<jetsize+1;m++) { //+1 for centre of jet (n steps is n+1 rings!)
            if((abs(i-k) == m && abs(j-l) <= m) || (abs(i-k) <= m && abs(j-l) == m)) { 
              //i.e. we are now in ring m
              if (m <= vetowindowsize) localsums[m] += ttArray[k][newl]; 
              areas[m]+=1;
              if(m == jetsize) { //i.e. we are in the outer ring and want to parameterise PU
                if( (k-i) == m && abs(j-l) <= (m-1) ) { outerstrips[0].first += ttArray[k][newl];outerstrips[0].second+=1;}
                if( (i-k) == m && abs(j-l) <= (m-1) ) { outerstrips[1].first += ttArray[k][newl];outerstrips[1].second+=1;}
                if( (l-j) == m && abs(i-k) <= (m-1) ) { outerstrips[2].first += ttArray[k][newl];outerstrips[2].second+=1;}
                if( (j-l) == m && abs(i-k) <= (m-1) ) { outerstrips[3].first += ttArray[k][newl];outerstrips[3].second+=1;}
              }

              if(m > 0 && m <= vetowindowsize) { //i.e. don't compare the central tower or towers outside vetowindowsize
                jetarea++;
                if((k+l) > (i+j) ) { if(ttArray[k][newl] > ttArray[i][j]) { numtowersaboveme++; } }
                else if( ((k+l) == (i+j)) && (k-i) > (l-j)) { if(ttArray[k][newl] > ttArray[i][j]) { numtowersaboveme++; } } //this line is to break the degeneracy along the diagonal treating top left different to bottom right
                else { if(ttArray[k][newl] >= ttArray[i][j]) { numtowersaboveme++; } }

              }
              break; //no point continuining since can only be a member of one ring
            }
          }

        }
      }

      //now we have a jet candidate centred at i,j, with the ring energies and areas defined

      //now we have the L1 jet candidate:
      if(numtowersaboveme == 0 )
      {
        //std::cout << ttArray[i][j] << "  "  << localsums[0] << std::endl;
        {
          if(ttArray[i][j] > seedthresh1) 
          {

            double totalenergy=0.0;
            //std::cout << "iEta: " << g.old_iEta(i) << ", iPhi: " << g.old_iPhi(j) << ", r0: " << localsums[0] <<  ", r1: " << localsums[1] << ", r2: " << localsums[2] << ", r3: " << localsums[3] << ", r4: " << localsums[4] << std::endl;
            for(int ring=0; ring < (int)localsums.size(); ring++) { totalenergy += localsums[ring]; }
            //this is with PUS:
            //for(int ring=0; ring < (int)localsums.size()-1; ring++) { totalenergy += localsums[ring]; }
            //std::sort(outerstrips.begin(),outerstrips.end());
            //totalenergy = totalenergy - (3.5 * (outerstrips[1] + outerstrips[2]));

            //this means we have a viable candidate
            if(totalenergy > 0.0) {
              //Store as a leaf candidate with charge of 0 and 0 mass
              //Multiply the energy by 0.5 to convert to GeV
              math::PtEtaPhiMLorentzVector tempJet(0.5*totalenergy, g.eta(g.old_iEta(i)), g.phi(g.old_iPhi(j)),0.);
              uncalibL1Jets.push_back(L1JetParticle(tempJet, L1JetParticle::JetType::kCentral,0));
              jetareas.push_back(jetarea);
              //THEY ARENT THE SAME
              //std::cout << "Eta: " << g.eta(g.old_iEta(i)) << "  Phi: " << g.phi(g.old_iPhi(j)) << std::endl;
              //std::cout << "Eta: " << uncalibL1Jets.back().eta() << "  Phi: " << uncalibL1Jets.back().phi() << std::endl;
            }
          }

        } 
      }
    }
  }

  //BUG HERE********************************
  //Perform global rho subtraction
  double median_energy = getMedian(uncalibL1Jets, jetareas);

  for(unsigned i =0; i<uncalibL1Jets.size(); i++){
    L1JetParticle newJet= L1JetParticle(math::PtEtaPhiMLorentzVector(uncalibL1Jets.at(i).pt()-median_energy*jetareas[i],uncalibL1Jets.at(i).eta(),uncalibL1Jets.at(i).phi(),0.), L1JetParticle::JetType::kCentral,0);
    
    //std::cout << "Eta: " << uncalibL1Jets.at(i).eta() << "  Phi: " << uncalibL1Jets.at(i).phi() << std::endl;

    uncalibL1Jets.at(i)=newJet;

    //THEYRE THE SAME NOW!
    //std::cout << "Eta: " << uncalibL1Jets.at(i).eta() << "  Phi: " << uncalibL1Jets.at(i).phi() << std::endl;

  }
  // *****************************************
  //sort by highest pT before ending
  std::sort(uncalibL1Jets.begin(), uncalibL1Jets.end(), sortbypt);  

  //Produce the calibrated l1 Jets, only calibrating down to 30 GeV now
  //std::auto_ptr<std::vector<reco::LeafCandidate> > l1Jets(new std::vector<reco::LeafCandidate>());
  std::vector<L1JetParticle> l1Jets;

  //Calibrations only work down to 30GeV as of now
  l1Jets = Stage2Calibrations::calibrateL1Jets(uncalibL1Jets,30.,9999.);

  //The thresholds are at least 30GeV due to calibrations
  //std::auto_ptr<std::vector<reco::LeafCandidate> > mht(new std::vector<reco::LeafCandidate>());
  //std::auto_ptr<double> ht(new double);
  std::vector<reco::LeafCandidate>  mht;
  mht.push_back(calculateMHT(l1Jets,mhtThreshold_));
  double ht = calculateHT(l1Jets,htThreshold_);

  std::auto_ptr<std::vector<reco::LeafCandidate> > mhtPtr( new std::vector<reco::LeafCandidate>() );
  std::auto_ptr<std::vector<L1JetParticle> > l1JetsPtr( new std::vector<L1JetParticle>() );
  std::auto_ptr<std::vector<L1JetParticle> > uncalibL1JetsPtr( new std::vector<L1JetParticle>() );
  std::auto_ptr<double> htPtr( new double() );

  *mhtPtr=mht;
  *l1JetsPtr=l1Jets;
  *uncalibL1JetsPtr=uncalibL1Jets;
  *htPtr=ht;

  iEvent.put(uncalibL1JetsPtr,"l1Stage2JetsUncalib");
  iEvent.put(l1JetsPtr,"l1Stage2Jets");
  iEvent.put(mhtPtr,"l1Stage2Mht");
  iEvent.put(htPtr,"l1Stage2Ht");

}

// ------------ method called once each job just before starting event loop  ------------
void Stage2JetProducer::beginJob() {
}

// ------------ method called once each job just after ending the event loop  ------------
void Stage2JetProducer::endJob() {
}

// ------------ method called when starting to processes a run  ------------
void Stage2JetProducer::beginRun(edm::Run&, edm::EventSetup const&) {
}

// ------------ method called when ending the processing of a run  ------------
void Stage2JetProducer::endRun(edm::Run&, edm::EventSetup const&) {
}

// ------------ method called when starting to processes a luminosity block  ------------
void Stage2JetProducer::beginLuminosityBlock(edm::LuminosityBlock&, 
    edm::EventSetup const&) {
}

// ------------ method called when ending the processing of a luminosity block  ------------
void Stage2JetProducer::endLuminosityBlock(edm::LuminosityBlock&, 
    edm::EventSetup const&) {
}

double Stage2JetProducer::getMedian(const std::vector<L1JetParticle>& jets, const std::vector<int>& areas)
{
  //std::sort(jets.begin(),jets.end(),sortbyrho);
  std::vector<L1JetParticle> jetSort=jets;

  //Scale the pt of all the jets by there areas
  for(unsigned i =0; i< jetSort.size(); i++){
    L1JetParticle newJet= L1JetParticle(math::PtEtaPhiMLorentzVector(jetSort[i].pt()/areas[i],jetSort[i].eta(),jetSort[i].phi(),0.), l1extra::L1JetParticle::JetType::kCentral, 0);
    jetSort[i]=newJet;
  }

  //sort by the rho
  std::sort(jetSort.begin(),jetSort.end(), sortbypt);  
  double median_energy=0.;
  if(jetSort.size() > 2) {
    if(jetSort.size() % 2 == 0 ) {
      int index1 = jetSort.size() / 2;
      median_energy=((double)((double)jetSort.at(index1).pt() + (double)jetSort.at(index1+1).pt()))/2.0;

    } else {
      int index1 = (int)jetSort.size() / 2;
      median_energy=(double)jetSort.at(index1).pt();
    }
  }
  return median_energy;

}

double Stage2JetProducer::calculateHT(const std::vector<L1JetParticle>& jets, const double& thresh) {
  double ht=0.0;
  for(unsigned int i=0; i< jets.size(); i++) {
    if (jets[i].pt() > thresh)  ht += jets[i].pt();
  }
  return ht;
}

LeafCandidate Stage2JetProducer::calculateMHT(const std::vector<L1JetParticle> & jets, const double& thresh) {

  double mht_x=0.0;
  double mht_y=0.0;
  for(unsigned int i=0; i< jets.size(); i++) {
    if (jets[i].pt() > thresh)
    {
      mht_x -= cos(jets[i].phi())*jets[i].pt();
      mht_y -= sin(jets[i].phi())*jets[i].pt();
    }
  }

  //double phi = atan2(mht_y,mht_x);
  LeafCandidate mht = LeafCandidate(0,math::XYZTLorentzVector(mht_x,mht_y,0.,sqrt(mht_x*mht_x+mht_y*mht_y)));
  //LeafCandidate mht = LeafCandidate(0,math::PtEtaPhiMLorentzVector(sqrt(mht_x*mht_x+mht_y*mht_y),0.,phi,0.), l1extra::L1JetParticle::JetType::kCentral, 0);

  return mht;
}



// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void Stage2JetProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
DEFINE_FWK_MODULE(Stage2JetProducer);

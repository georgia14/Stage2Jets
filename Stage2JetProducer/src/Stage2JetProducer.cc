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
  produces<std::vector<L1JetParticle> >("l1Stage2JetsNoPUS");
  produces<std::vector<L1JetParticle> >("l1Stage2JetsNoPUSUncalib");
  produces<std::vector<L1JetParticle> >("l1Stage2JetsDonutPUS");
  produces<std::vector<L1JetParticle> >("l1Stage2JetsDonutPUSUncalib");
  produces<std::vector<L1JetParticle> >("l1Stage2JetsGlobalPUS");
  produces<std::vector<L1JetParticle> >("l1Stage2JetsGlobalPUSUncalib");
  produces<std::vector<L1EtMissParticle> >("l1Stage2NoPUSMht");
  produces<std::vector<L1EtMissParticle> >("l1Stage2DonutPUSMht");
  produces<std::vector<L1EtMissParticle> >("l1Stage2GlobalPUSMht");
  produces<std::vector<L1EtMissParticle> >("l1Stage2Met");
  //produces<double>("l1Stage2Ht");

  //Get the configs
  towersTag_ = iConfig.getParameter<edm::InputTag>("towerInput");
  mhtThreshold_ = iConfig.getParameter<double>("mhtThreshold");
  htThreshold_ = iConfig.getParameter<double>("htThreshold");
  jetEtaCut_ = iConfig.getParameter<double>("jetEtaCut");

}


Stage2JetProducer::~Stage2JetProducer() {
}

// ------------ method called to produce the data  ------------
void Stage2JetProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  //Get the trigger towers and put them in an array
  edm::Handle<l1slhc::L1CaloTowerCollection> triggerTowers;
  iEvent.getByLabel(towersTag_, triggerTowers);

  //Maybe make these configurable but not for now
  //int jetsize=5;
  //int vetowindowsize=4;
  int seedthresh1=0;
  int seedthresh2=0;
  int nstrips=4;


  //Make these incase you want them in the future
  double ET=0.;
  double met_x=0.;
  double met_y=0.;

  std::vector< std::vector<int> > input(56, std::vector<int>(72, 0)); //this is just a container for the (E+H) per tower

  //std::auto_ptr<std::vector<reco::LeafCandidate> > uncalibL1Jets(new
  //    std::vector<reco::LeafCandidate>());
  std::vector<L1JetParticle> uncalibL1Jets;
  std::vector<L1JetParticle> uncalibChunkyJets;

  for(auto j=triggerTowers->begin(); j!=triggerTowers->end(); j++) {

    if ( abs((*j).iEta()) > 28 ) { continue; } //i.e. |eta| < 3 only

    input[g.new_iEta((*j).iEta())][g.new_iPhi((*j).iPhi())] = ((*j).E() + (*j).H());
    ET += (*j).E()+(*j).H();
    met_x -= cos(g.phi((*j).iPhi())) * ((*j).E() + (*j).H());
    met_y -= sin(g.phi((*j).iPhi())) * ((*j).E() + (*j).H());
    //so now ttArray is on the scale ieta 0-56, iphi 0-71
  }

  //Put in Mask for chunky donut
  std::vector<std::vector<int> > mask, mask_donut;

  //-------------The Masks------------------//

  mask=mask_square_9by9();
  mask_donut= mask_chunky_donut_15by15();

  //----------------------------------------//

  std::vector<double> jetareas; //to hold the ring areas (i.e. when we get up against the boundaries)
  std::vector<int> seedValue; //to hold the ring areas (i.e. when we get up against the boundaries)

  TriggerTowerGeometry g;
  int etasize=(mask.size()-1)/2;
  int phisize=(mask.at(0).size()-1)/2;
  int nringsveto =etasize;

  int etasizedonut=(mask_donut.size()-1)/2;
  int phisizedonut=(mask_donut.at(0).size()-1)/2;
  int nstripsdonut = nstrips;
  /////////////////
  //std::cout << input.size() << ", " << input[0].size() << std::endl;
  for ( int i = 0; i < (int)input.size(); i++) {
    for ( int j = 0; j < (int)input[i].size(); j++) {
      std::vector<uint8_t> jetTower(81,0);
      //int seedtower = input[i][j];  
      //std::cout << "new: (" << i << ", " << j << ", " << input[i][j] << ")" << std::endl;
      int numtowersaboveme=0;
      int numtowersabovethresh=0;

      std::vector<int> localsums(nringsveto+1,0); //to hold the ring sums (+1 for centre)
      std::vector<int> localmax(nringsveto+1,0); //to hold the ring sums (+1 for centre)
      std::vector<double> areas(nringsveto+1,0.); //to hold the ring areas (i.e. when we get up against the boundaries)
      //std::vector<int> outerstrips(nstripsdonut,0); //to hold the energies in the 4 surrounding outer strips (excluding corners)
      std::vector<std::pair<int,double>> outerstrips(nstripsdonut,std::make_pair(0,0.)); //to hold the energies in the 4 surrounding outer strips (excluding corners)//AND THEIR AREAS!
      //areas[0]=1;
      areas[0]=g.towerEtaSize(g.old_iEta(i));;
      //double jetarea = 1;
      double jetarea = g.towerEtaSize(g.old_iEta(i));;
      //int pusarea=0;
      for(int l=(j-phisizedonut); l<=(j+phisizedonut); l++) {
        for(int k=(i-etasizedonut); k<=(i+etasizedonut); k++) {
          if(k < 0 || k > 55) continue; //i.e. out of bounds of eta<3
          //std::cout << " k = " << k << ", l = " << l << ", i =" << i << ", j = " << j << std::endl;
          unsigned int dk = k-i+etasizedonut;
          unsigned int dl = l-j+phisizedonut;


          //	       std::cout << dk << std::endl;
          //	       std::cout << dl << std::endl;
          //make a co-ordinate transform at the phi boundary
          int newl;
          if(l < 0) { newl = l+72; } 
          else if (l > 71) { newl = l-72; } 
          else { newl = l; }

          //if (l != j && k != i)
          if(input[k][newl] > seedthresh2) { numtowersabovethresh++; }
          if(dl < mask_donut.size() && dk < mask_donut.at(0).size())
          {
            if(mask_donut[dl][dk] != 0)
            {
              double towerArea = g.towerEtaSize(g.old_iEta(k));
              outerstrips[mask_donut[dl][dk]-1].first+=input[k][newl];
              outerstrips[mask_donut[dl][dk]-1].second+=towerArea;
              //	     std::cout << mask_donut[dl][dk];
            }
          }

        }
        //std::cout << std::endl;
      }
      //std::cout << outerstrips[0].first <<" "<< outerstrips[0].second << std::endl;
      //std::cout << std::endl;
      double jetSecMomEta = 0.; 
      double jetSecMomPhi = 0.; 
      double jetFirMomEta = 0.; 
      double jetFirMomPhi = 0.; 
      double jetCovEtaPhi = 0.;
      int towIndex = 0;
      for(int l=(j-phisize); l<=(j+phisize); l++) {
        for(int k=(i-etasize); k<=(i+etasize); k++) {
          unsigned int dk = k-i+etasize;
          unsigned int dl = l-j+phisize;
          if(k < 0 || k > 55) continue; //i.e. out of bounds of eta<3
          //std::cout << " k = " << k << ", l = " << l << ", i =" << i << ", j = " << j << std::endl;
          //	       std::cout << dk << std::endl;
          //	       std::cout << dl << std::endl;
          //make a co-ordinate transform at the phi boundary
          int newl;
          if(l < 0) { newl = l+72; } 
          else if (l > 71) { newl = l-72; } 
          else { newl = l; }

          //if (l != j && k != i)
          {
            jetTower.at(towIndex)=(input[k][newl]);
          }
          towIndex++;
          if(input[k][newl] > seedthresh2) { numtowersabovethresh++; }
          if(dl < mask.size() && dk < mask.at(0).size())
          {
            jetSecMomPhi+=((int)dl-phisize)*((int)dl-phisize)*input[k][newl]*input[k][newl];
            jetSecMomEta+=((int)dk-etasize)*((int)dk-etasize)*input[k][newl]*input[k][newl];
            jetFirMomPhi+=((int)dl-phisize)*input[k][newl];
            jetFirMomEta+=((int)dk-etasize)*input[k][newl];
            jetCovEtaPhi+=((int)dk-etasize)*((int)dl-phisize)*input[k][newl]*input[k][newl];
            if (mask[dl][dk] == 2){if(input[k][newl]>input[i][j]) {numtowersaboveme++;}}
            else if (mask[dl][dk] == 1){if(input[k][newl]>=input[i][j]) {numtowersaboveme++;}}
            /*
               if((k+l) > (i+j) ) { if(input[k][newl] > input[i][j]) { numtowersaboveme++; } }
               else if( ((k+l) == (i+j)) && (k-i) > (l-j)) { if(input[k][newl] > input[i][j]) { numtowersaboveme++; } } //this line is to break the degeneracy along the diagonal treating top left different to bottom right
               else { if(input[k][newl] >= input[i][j]) { numtowersaboveme++; } }
               */
            for( int m=0; m<nringsveto+1;m++) { //+1 for centre of jet (n steps is n+1 rings!)
              if((abs(i-k) == m && abs(j-l) <= m) || (abs(i-k) <= m && abs(j-l) == m)) { 
                //i.e. we are now in ring m
                if(input[k][newl]>localmax[m]) localmax[m] = input[k][newl];
                localsums[m] += input[k][newl]; 
                double towerArea = g.towerEtaSize(g.old_iEta(k));
                if (mask[dl][dk] != 0) {areas[m] += towerArea; jetarea+=towerArea;}
                break; //no point continuining since can only be a member of one ring
              }
            }
          }

          if (numtowersaboveme > 0) break;
        }
        if (numtowersaboveme > 0) break;
      }
      //now we have a jet candidate centred at i,j, with the ring energies and areas defined
      //now we have the L1 jet candidate:
      if(numtowersaboveme == 0 && input[i][j] > seedthresh1) {
        double totalenergy=0.0;
        //std::cout << "iEta: " << g.old_iEta(i) << ", iPhi: " << g.old_iPhi(j) << ", r0: " << localsums[0] <<  ", r1: " << localsums[1] << ", r2: " << localsums[2] << ", r3: " << localsums[3] << ", r4: " << localsums[4] << std::endl;
        for(int ring=0; ring < (int)localsums.size(); ring++) { totalenergy += localsums[ring]; }
        jetSecMomEta = jetSecMomEta/(totalenergy*totalenergy);
        jetSecMomPhi = jetSecMomPhi/(totalenergy*totalenergy);
        jetFirMomEta = jetFirMomEta/(totalenergy);
        jetFirMomPhi = jetFirMomPhi/(totalenergy);
        //this is with PUS:
        if(totalenergy > 0.0 && fabs(g.eta(g.old_iEta(i)))<jetEtaCut_  ) {
          //L1_jJets.push_back(jJet(totalenergy, g.old_iEta(i), g.old_iPhi(j), localsums,localmax,jetFirMomEta,jetFirMomPhi,jetSecMomEta,jetSecMomPhi,jetCovEtaPhi, areas, outerstrips,jetTower,jetarea));

          math::PtEtaPhiMLorentzVector tempJet(0.5*totalenergy, g.eta(g.old_iEta(i)), g.phi(g.old_iPhi(j)),0.);
          uncalibL1Jets.push_back(L1JetParticle(tempJet, L1JetParticle::JetType::kCentral,0));
          jetareas.push_back(jetarea);
          seedValue.push_back(localsums[0]);

          //Apply seed threshold
          if(localsums[0]>5){

            //Do chunky subtraction
            std::sort(outerstrips.begin(),outerstrips.end());
            double chunkyEnergy = totalenergy- 
              jetarea*(outerstrips[1].first+outerstrips[2].first)/(outerstrips[1].second+outerstrips[2].second);

            if(chunkyEnergy>0.0){
              math::PtEtaPhiMLorentzVector tempJet2(0.5*chunkyEnergy, g.eta(g.old_iEta(i)), g.phi(g.old_iPhi(j)),0.);
              uncalibChunkyJets.push_back(L1JetParticle(tempJet2, L1JetParticle::JetType::kCentral,0));

            }

          }

        }

      }

    }
  }

  //
  //---------- Do pileup subtraction jets -----------------------------//

  //Do global rho subtraction
  //
  std::vector<L1JetParticle> uncalibGlobalJets;

  double median_energy = getMedian(uncalibL1Jets, jetareas);

  for(unsigned i =0; i<uncalibL1Jets.size(); i++){
    L1JetParticle newJet= L1JetParticle(math::PtEtaPhiMLorentzVector(uncalibL1Jets.at(i).pt()-median_energy*jetareas[i],uncalibL1Jets.at(i).eta(),uncalibL1Jets.at(i).phi(),0.), L1JetParticle::JetType::kCentral,0);

    //Check the seed
    if(seedValue[i]>5 && uncalibL1Jets.at(i).pt()-median_energy > 0.){
      uncalibGlobalJets.push_back(newJet);
    }
  }

  //Sort the global jets
  std::sort(uncalibGlobalJets.begin(),uncalibGlobalJets.end(),sortbypt);

  //Calibrate the global jets
  std::vector<L1JetParticle> globalJets = Stage2Calibrations::calibrateL1Jets(uncalibGlobalJets,"global",10.,9999.);
  std::sort(globalJets.begin(),globalJets.end(),sortbypt);

  //Put into the event
  std::auto_ptr<std::vector<L1JetParticle> > globalJetsPtr( new std::vector<L1JetParticle>() );
  std::auto_ptr<std::vector<L1JetParticle> > uncalibGlobalJetsPtr( new std::vector<L1JetParticle>() );
  *globalJetsPtr=globalJets;
  *uncalibGlobalJetsPtr=uncalibGlobalJets;

  iEvent.put(uncalibGlobalJetsPtr,"l1Stage2JetsGlobalPUSUncalib");
  iEvent.put(globalJetsPtr,"l1Stage2JetsGlobalPUS");

  //Do chunky donut subtraction
  //
  //Sort the chunky jets
  std::sort(uncalibChunkyJets.begin(),uncalibChunkyJets.end(),sortbypt);

  //Calibrate the chunky jets
  std::vector<L1JetParticle> chunkyJets = Stage2Calibrations::calibrateL1Jets(uncalibChunkyJets,"chunky",10.,9999.);

  //Sort the calibrated jets
  std::sort(chunkyJets.begin(),chunkyJets.end(),sortbypt);

  //Put into the event
  std::auto_ptr<std::vector<L1JetParticle> > chunkyJetsPtr( new std::vector<L1JetParticle>() );
  std::auto_ptr<std::vector<L1JetParticle> > uncalibChunkyJetsPtr( new std::vector<L1JetParticle>() );
  *chunkyJetsPtr=chunkyJets;
  *uncalibChunkyJetsPtr=uncalibChunkyJets;

  iEvent.put(uncalibChunkyJetsPtr,"l1Stage2JetsDonutPUSUncalib");
  iEvent.put(chunkyJetsPtr,"l1Stage2JetsDonutPUS");

  //----------------------------------------------------------------//
  //
  //---------- Do for unsubtracted jets -----------------------------//

  //sort by highest pT before ending
  std::sort(uncalibL1Jets.begin(), uncalibL1Jets.end(), sortbypt);  

  //Produce the calibrated l1 Jets, only calibrating down to 30 GeV now
  std::vector<L1JetParticle> l1Jets;

  //Calibrations only work down to 30GeV as of now
  l1Jets = Stage2Calibrations::calibrateL1Jets(uncalibL1Jets,"nopus",10.,9999.);

  //Sort the calibrated jets
  std::sort(l1Jets.begin(), l1Jets.end(), sortbypt);  

  //Put into the event
  std::auto_ptr<std::vector<L1JetParticle> > l1JetsPtr( new std::vector<L1JetParticle>() );
  std::auto_ptr<std::vector<L1JetParticle> > uncalibL1JetsPtr( new std::vector<L1JetParticle>() );

  *l1JetsPtr=l1Jets;
  *uncalibL1JetsPtr=uncalibL1Jets;

  iEvent.put(uncalibL1JetsPtr,"l1Stage2JetsNoPUSUncalib");
  iEvent.put(l1JetsPtr,"l1Stage2JetsNoPUS");

  //----------------------------------------------------------------//
  //
  //---------------------Make the Energy sums-------------------------//

  std::auto_ptr<std::vector<L1EtMissParticle> > mhtPtr( new std::vector<L1EtMissParticle>() );
  std::auto_ptr<std::vector<L1EtMissParticle> > metPtr( new std::vector<L1EtMissParticle>() );
  std::auto_ptr<std::vector<L1EtMissParticle> > mhtChunkyPtr( new std::vector<L1EtMissParticle>() );
  std::auto_ptr<std::vector<L1EtMissParticle> > mhtGlobalPtr( new std::vector<L1EtMissParticle>() );
  std::auto_ptr<double> htPtr( new double() );

  //Vanilla jets
  std::vector<L1EtMissParticle>  mht;
  mht.push_back(calculateMHT(l1Jets,mhtThreshold_,htThreshold_));
  //double ht = calculateHT(l1Jets,htThreshold_);

  std::vector<L1EtMissParticle>  mhtChunky;
  mhtChunky.push_back(calculateMHT(chunkyJets,mhtThreshold_,htThreshold_));

  std::vector<L1EtMissParticle>  mhtGlobal;
  mhtGlobal.push_back(calculateMHT(globalJets,mhtThreshold_,htThreshold_));

  std::vector<L1EtMissParticle> met;
  met.push_back(
      L1EtMissParticle(math::PtEtaPhiMLorentzVector(sqrt(met_x*met_x+met_y*met_y),0.,atan2(met_y,met_x),0.),
        L1EtMissParticle::EtMissType::kMET,ET));

  *mhtPtr=mht;
  *mhtChunkyPtr=mhtChunky;
  *mhtGlobalPtr=mhtGlobal;
  *metPtr=met;
  //*htPtr=ht;

  iEvent.put(mhtPtr,"l1Stage2NoPUSMht");
  iEvent.put(metPtr,"l1Stage2Met");
  iEvent.put(mhtChunkyPtr,"l1Stage2DonutPUSMht");
  iEvent.put(mhtGlobalPtr,"l1Stage2GlobalPUSMht");

  //----------------------------------------------------------------//
  //
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

double Stage2JetProducer::getMedian(const std::vector<L1JetParticle>& jets, const std::vector<double>& areas)
{
  //std::sort(jets.begin(),jets.end(),sortbyrho);
  std::vector<L1JetParticle> jetSort=jets;

  //Scale the pt of all the jets by their areas
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

L1EtMissParticle Stage2JetProducer::calculateMHT(const std::vector<L1JetParticle> & jets, const double& mhtThresh, const double& htThresh){

  double mht_x=0.0;
  double mht_y=0.0;
  double ht=0.0;
  for(unsigned int i=0; i< jets.size(); i++) {

    if (jets[i].pt() > htThresh)  ht += jets[i].pt();
    if (jets[i].pt() > mhtThresh)
    {
      mht_x -= cos(jets[i].phi())*jets[i].pt();
      mht_y -= sin(jets[i].phi())*jets[i].pt();
    }
  }

  double phi = atan2(mht_y,mht_x);
  //LeafCandidate mht = LeafCandidate(0,math::XYZTLorentzVector(mht_x,mht_y,0.,sqrt(mht_x*mht_x+mht_y*mht_y)));
  //LeafCandidate mht = LeafCandidate(0,math::PtEtaPhiMLorentzVector(sqrt(mht_x*mht_x+mht_y*mht_y),0.,phi,0.), l1extra::L1JetParticle::JetType::kCentral, 0);
  L1EtMissParticle mht = L1EtMissParticle(math::PtEtaPhiMLorentzVector(sqrt(mht_x*mht_x+mht_y*mht_y),0.,phi,0.),L1EtMissParticle::EtMissType::kMHT,ht);

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

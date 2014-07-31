#ifndef STAGE2JETPRODUCER_H
#define STAGE2JETPRODUCER_H

#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Candidate/interface/LeafCandidate.h"
#include "DataFormats/L1Trigger/interface/L1JetParticle.h"
#include "DataFormats/L1Trigger/interface/L1EtMissParticle.h"
#include "DataFormats/CaloTowers/interface/CaloTower.h"
#include "DataFormats/CaloTowers/interface/CaloTowerFwd.h"
#include "SimDataFormats/SLHC/interface/L1CaloTower.h"
#include "SimDataFormats/SLHC/interface/L1CaloTowerFwd.h"

#include "Stage2Jets/Stage2JetProducer/plugins/TriggerTowerGeometry_new.cc"
#include "Stage2Jets/Stage2JetProducer/plugins/CalibrationFunctions.cc"
#include "Stage2Jets/Stage2JetProducer/interface/mask.hh"

using namespace reco;
using namespace l1extra;

class Stage2JetProducer : public edm::EDProducer {
   public:
      explicit Stage2JetProducer(const edm::ParameterSet&);
      ~Stage2JetProducer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
      
      double getMedian(const std::vector<L1JetParticle>& jets, const std::vector<double>& areas);
      double calculateHT(const std::vector<L1JetParticle>& jets, const double& threshold);
      L1EtMissParticle calculateMHT(const std::vector<L1JetParticle>& jets, const double& threshold);

   private:
      virtual void beginJob() ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      virtual void beginRun(edm::Run&, edm::EventSetup const&);
      virtual void endRun(edm::Run&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);



      TriggerTowerGeometry g;

      //Define the input tags
      edm::InputTag towersTag_;

      double mhtThreshold_;
      double htThreshold_;
      double jetEtaCut_;

};

#endif

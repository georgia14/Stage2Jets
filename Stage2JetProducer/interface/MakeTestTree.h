#ifndef MAKETESTTREE_H
#define MAKETESTTREE_H

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "TTree.h"

#include "DataFormats/Candidate/interface/LeafCandidate.h"
#include "DataFormats/L1Trigger/interface/L1JetParticle.h"
#include "DataFormats/L1Trigger/interface/L1EtMissParticle.h"

#include "DataFormats/L1TCalorimeter/interface/CaloTower.h"
#include "DataFormats/L1TCalorimeter/interface/CaloCluster.h"

#include "DataFormats/CaloTowers/interface/CaloTower.h"
#include "DataFormats/CaloTowers/interface/CaloTowerFwd.h"
#include "SimDataFormats/SLHC/interface/L1CaloTower.h"
#include "SimDataFormats/SLHC/interface/L1CaloTowerFwd.h"

#include "Stage2Jets/Stage2JetProducer/plugins/TriggerTowerGeometry_new.cc"

using namespace l1t;
using namespace l1extra;

class MakeTestTree : public edm::EDAnalyzer {
   public:
      explicit MakeTestTree(const edm::ParameterSet&);
      ~MakeTestTree();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

      
   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

      TTree * tree;
      unsigned run;
      unsigned event;

      int Njet;
      int NjetUncalib;
      int NjetGlobal;
      int NjetGlobalUncalib;
      int NjetNoPus;
      int NjetNoPusUncalib;

      std::vector<int> *towerPt;
      std::vector<int> *towerEta;
      std::vector<int> *towerPhi;
      std::vector<int> *towerEtEm;
      std::vector<int> *towerEtHad;

      std::vector<double> *jetPt;
      std::vector<double> *jetPhi;
      std::vector<double> *jetEta;
      std::vector<double> *jetPhi_ph;
      std::vector<double> *jetEta_ph;
      std::vector<double> *jetUncalibPt;
      std::vector<double> *jetUncalibPhi;
      std::vector<double> *jetUncalibEta;
      std::vector<double> *jetUncalibPhi_ph;
      std::vector<double> *jetUncalibEta_ph;
      std::vector<double> *jetGlobalPt;
      std::vector<double> *jetGlobalPhi;
      std::vector<double> *jetGlobalEta;
      std::vector<double> *jetGlobalPhi_ph;
      std::vector<double> *jetGlobalEta_ph;
      std::vector<double> *jetGlobalUncalibPt;
      std::vector<double> *jetGlobalUncalibPhi;
      std::vector<double> *jetGlobalUncalibEta;
      std::vector<double> *jetGlobalUncalibPhi_ph;
      std::vector<double> *jetGlobalUncalibEta_ph;
      std::vector<double> *jetNoPusPt;
      std::vector<double> *jetNoPusPhi;
      std::vector<double> *jetNoPusEta;
      std::vector<double> *jetNoPusPhi_ph;
      std::vector<double> *jetNoPusEta_ph;
      std::vector<double> *jetNoPusUncalibPt;
      std::vector<double> *jetNoPusUncalibPhi;
      std::vector<double> *jetNoPusUncalibEta;
      std::vector<double> *jetNoPusUncalibPhi_ph;
      std::vector<double> *jetNoPusUncalibEta_ph;

      double mht;
      double mhtPhi;
      double ht;
      double met;
      double metPhi;
      double et;
      double mhtGlobal;
      double mhtPhiGlobal;
      double htGlobal;
      double mhtNoPus;
      double mhtPhiNoPus;
      double htNoPus;

      TriggerTowerGeometry g;

      edm::InputTag towersTag_;
      edm::EDGetToken m_towerToken;

};
#endif

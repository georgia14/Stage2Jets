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
      std::vector<double> *jetPt;
      std::vector<double> *jetPhi;
      std::vector<double> *jetEta;
      std::vector<double> *jetUncalibPt;
      std::vector<double> *jetUncalibPhi;
      std::vector<double> *jetUncalibEta;
      std::vector<double> *jetGlobalPt;
      std::vector<double> *jetGlobalPhi;
      std::vector<double> *jetGlobalEta;
      std::vector<double> *jetGlobalUncalibPt;
      std::vector<double> *jetGlobalUncalibPhi;
      std::vector<double> *jetGlobalUncalibEta;
      std::vector<double> *jetNoPusPt;
      std::vector<double> *jetNoPusPhi;
      std::vector<double> *jetNoPusEta;
      std::vector<double> *jetNoPusUncalibPt;
      std::vector<double> *jetNoPusUncalibPhi;
      std::vector<double> *jetNoPusUncalibEta;
      double mht;
      double mhtPhi;
      double ht;
      double mhtGlobal;
      double mhtPhiGlobal;
      double htGlobal;
      double mhtNoPus;
      double mhtPhiNoPus;
      double htNoPus;

};
#endif

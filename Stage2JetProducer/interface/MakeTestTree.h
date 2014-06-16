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
      double mht;
      double mhtPhi;
      double ht;

};
#endif

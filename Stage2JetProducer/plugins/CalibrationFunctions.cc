#ifndef CALIBRATIONFUNCTIONS_HH
#define CALIBRATIONFUNCTIONS_HH

#include "DataFormats/Candidate/interface/LeafCandidate.h"
#include "Stage2Jets/Stage2JetProducer/plugins/TriggerTowerGeometry_new.cc"
#include "DataFormats/L1Trigger/interface/L1JetParticle.h"

using namespace reco;
using namespace l1extra;

namespace Stage2Calibrations{

  double calibFit( Double_t *v, Double_t *par ){

    // JETMET uses log10 rather than the ln used here...
    double logX = log(v[0]);

    double term1 = par[1] / ( logX * logX + par[2] );
    double term2 = par[3] * exp( -par[4]*((logX - par[5])*(logX - par[5])) );

    // Final fitting function 
    double f    = par[0] + term1 + term2; 

    return f;
    //return 1.0/f;
  }


  std::vector<L1JetParticle> calibrateL1Jets(const std::vector<L1JetParticle>& inJets, double ptMin, double ptMax){

    //Look up table for global PUS
    double lut[48]= //Goes down to 8
    {
      0.906907,       9.799444,       4.898061,       1.048949,       0.136667,       0.254435,
      0.979404,       7.034494,       3.436158,       2.262959,       0.109428,       0.247115,
      0.923403,       9.999900,       4.586094,       1.717475,       0.104799,       0.000067,
      0.863852,       9.999617,       4.358355,       1.252304,       0.091954,       0.000096,
      0.738636,       9.968689,       1.752450,       0.769343,       0.050463,       0.004559,
      0.836676,       9.999992,       2.904923,       1.273807,       0.077622,       0.000002,
      0.845303,       9.985136,       2.998129,       1.810531,       0.084673,       0.002361,
      0.923338,       8.541081,       5.108460,       1.285704,       0.109812,       0.003085,
    };


    TriggerTowerGeometry g;

    std::vector<L1JetParticle> outJets;

    for(auto iJet = inJets.begin(); iJet!=inJets.end(); iJet++){

      //If the pt of the jet is outside of the calibration range, add the jet to the calibrated jets and do nothing
      if(iJet->pt()>ptMax){
        outJets.push_back(*iJet);
      }
      else if(iJet->pt()>ptMin){

        double p[6]; //These are the parameters of the fit
        double v[1]; //This is the pt value

        //Load the lut based on the correct eta bin
        if(g.iEta(iJet->eta())>=-28 && g.iEta(iJet->eta())<-21){
          p[0]=lut[0];
          p[1]=lut[1];
          p[2]=lut[2];
          p[3]=lut[3];
          p[4]=lut[4];
          p[5]=lut[5];
        }else if(g.iEta(iJet->eta())>=-21 && g.iEta(iJet->eta())<-14){
          p[0]=lut[6];
          p[1]=lut[7];
          p[2]=lut[8];
          p[3]=lut[9];
          p[4]=lut[10];
          p[5]=lut[11];
        }else if(g.iEta(iJet->eta())>=-14 && g.iEta(iJet->eta())<7){
          p[0]=lut[12];
          p[1]=lut[13];
          p[2]=lut[14];
          p[3]=lut[15];
          p[4]=lut[16];
          p[5]=lut[17];
        }else if(g.iEta(iJet->eta())>=7 && g.iEta(iJet->eta())<0){
          p[0]=lut[18];
          p[1]=lut[19];
          p[2]=lut[20];
          p[3]=lut[21];
          p[4]=lut[22];
          p[5]=lut[23];
        }else if(g.iEta(iJet->eta())>0 && g.iEta(iJet->eta())<=7){
          p[0]=lut[24];
          p[1]=lut[25];
          p[2]=lut[26];
          p[3]=lut[27];
          p[4]=lut[28];
          p[5]=lut[29];
        }else if(g.iEta(iJet->eta())>7 && g.iEta(iJet->eta())<=14){
          p[0]=lut[30];
          p[1]=lut[31];
          p[2]=lut[32];
          p[3]=lut[33];
          p[4]=lut[34];
          p[5]=lut[35];
        }else if(g.iEta(iJet->eta())>14 && g.iEta(iJet->eta())<=21){
          p[0]=lut[36];
          p[1]=lut[37];
          p[2]=lut[38];
          p[3]=lut[39];
          p[4]=lut[40];
          p[5]=lut[41];
        }else if(g.iEta(iJet->eta())>21 && g.iEta(iJet->eta())<=28){
          p[0]=lut[42];
          p[1]=lut[43];
          p[2]=lut[44];
          p[3]=lut[45];
          p[4]=lut[46];
          p[5]=lut[47];
        }
        v[0]=iJet->pt();
        double correction=1.0*calibFit(v,p);

        L1JetParticle newJet= L1JetParticle(math::PtEtaPhiMLorentzVector(correction*iJet->pt(),iJet->eta(),iJet->phi(),0.), l1extra::L1JetParticle::JetType::kCentral, 0);

        outJets.push_back(newJet);

      }

    }
    return outJets;
  }

}
#endif


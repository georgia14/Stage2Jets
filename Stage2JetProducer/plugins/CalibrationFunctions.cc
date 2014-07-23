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


  std::vector<L1JetParticle> calibrateL1Jets(const std::vector<L1JetParticle>& inJets,TString type, double ptMin, double ptMax){

    //Look up table for global PUS
    double lut[48];

    if(type=="nopus"){

      lut=
      {
        1.094573,	1.914835,	-0.999952,	0.529477,	0.198517,	2.045523,
        1.117407,	3.225482,	-0.971224,	0.772645,	0.218468,	2.237276,
        1.113721,	1.907837,	0.386244,	0.914753,	0.162511,	1.750241,
        1.021381,	5.259399,	1.337182,	0.351027,	0.276796,	2.730577,
        1.008247,	5.990503,	9.975187,	0.731619,	0.156460,	1.633716,
        1.133472,	1.367577,	-0.368074,	0.932853,	0.168657,	1.787382,
        1.082999,	5.082858,	-0.170875,	0.588396,	0.252840,	2.479732,
        1.090921,	2.443421,	-0.996960,	0.438731,	0.223473,	2.294779,
      };

    }else if(type=="chunky"){

      lut=
      {
        1.004215,	7.236275,	10.000000,	0.334575,	0.521609,	2.836302,
        1.159529,	1.514186,	6.953342,	1.515749,	0.201221,	1.685506,
        1.263089,	0.000000,	9.963021,	1.664812,	0.233152,	1.841630,
        1.149536,	0.981542,	-0.998390,	1.522643,	0.177708,	1.449095,
        1.164966,	0.099319,	-0.913029,	1.830581,	0.157528,	1.169878,
        1.260835,	0.000000,	9.923617,	1.685096,	0.225176,	1.792198,
        1.198797,	0.048253,	2.740067,	1.662546,	0.188187,	1.585532,
        1.085414,	4.114091,	10.000000,	0.439377,	0.432024,	2.703808,
      };

    }else{
      std::cout << "Invalid jet type for calibration, not calibrating" << std::endl;
      return std::vector<L1JetParticle>();
    }

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


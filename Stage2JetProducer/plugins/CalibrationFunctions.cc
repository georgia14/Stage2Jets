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

    double nopusLut[48]=
    {
      1.117222,	0.115154,	1.910615,	1.013722,	0.141106,	0.932544,
      1.097168,	3.065602,	-1.000000,	0.769128,	0.237118,	2.149065,
      1.136084,	4.306006,	-1.000000,	0.560253,	0.387929,	2.712581,
      1.120055,	1.720310,	-1.000000,	0.701081,	0.236731,	2.260849,
      1.147088,	0.000000,	2.021120,	1.163963,	0.157996,	1.465768,
      1.126069,	4.677360,	-1.000000,	0.507453,	0.411144,	2.794454,
      1.089151,	3.979116,	-1.000000,	0.595467,	0.288850,	2.478056,
      1.135149,	0.010747,	7.372626,	0.961446,	0.153867,	1.104507,
    };


    double chunkyLut[48]=
    {
      1.162001,	0.000000,	0.506532,	0.660991,	0.313898,	2.431404,
      1.195166,	0.000000,	3.462064,	1.500481,	0.217480,	1.860866,
      1.268926,	0.000000,	9.016650,	1.562618,	0.264561,	2.044610,
      1.162806,	0.703037,	-0.998253,	1.440402,	0.203554,	1.713803,
      1.168191,	0.339318,	-0.925540,	1.539332,	0.191717,	1.600152,
      1.267062,	0.000000,	-0.951295,	1.551194,	0.259202,	2.032236,
      1.205309,	0.000248,	7.230826,	1.517030,	0.214394,	1.842128,
      1.174576,	0.156892,	9.999552,	0.635096,	0.343184,	2.509776,
    };

    double globalLut[48]=
    {
      1.134861,	0.219994,	7.507336,	1.170671,	0.188135,	1.429918,
      1.177166,	0.037915,	5.775766,	1.900280,	0.173129,	1.340181,
      1.255825,	0.000000,	5.354664,	1.505212,	0.247060,	1.963272,
      1.163506,	0.000068,	3.125683,	1.431822,	0.178695,	1.541276,
      1.159526,	0.000000,	4.377770,	1.462514,	0.170831,	1.454168,
      1.253616,	0.000000,	0.343636,	1.495963,	0.243225,	1.953174,
      1.176683,	0.550961,	-0.935741,	1.721197,	0.183204,	1.484217,
      1.156831,	0.000000,	2.441199,	1.104090,	0.208981,	1.608509,
    };

    double* lut;
    if(type=="nopus"){
      lut=nopusLut;
    }else if(type=="chunky"){
      lut=chunkyLut;
    }else if(type=="global"){
      lut=globalLut;
    }else{
      std::cout << "Invalid jet type for L1 Jet calibration, not calibrating" << std::endl;
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


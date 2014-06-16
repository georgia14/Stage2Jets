#ifndef TRIGGER_TOWER_GEO_NEW
#define TRIGGER_TOWER_GEO_NEW

#include <TMath.h>

//modified JM, AR
//
//Used to get the correct mapping for the Trigger Towers

class TriggerTowerGeometry
{

  private:
    double mSizeEta[32];
    double mMappingEta[32];

    const double mPhiSize;

  public:
    TriggerTowerGeometry() :
      mPhiSize( TMath::TwoPi()/72.0 ) {

        //std::cout << "Trigger Tower geometry:" << std::endl; //for debug

        for ( int i = 0; i != 20; ++i ) {
          mSizeEta[i] = 0.087; //mSizeEta defines the divisions of each of the bins
        }

        mSizeEta[20] = 0.09;
        mSizeEta[21] = 0.1;
        mSizeEta[22] = 0.113;
        mSizeEta[23] = 0.129;
        mSizeEta[24] = 0.15;
        mSizeEta[25] = 0.178;
        mSizeEta[26] = 0.15;
        mSizeEta[27] = 0.35;

        for ( int i = 28; i != 32 ; ++i ) {
          mSizeEta[i] = 0.5; //the famous HF...
        }

        double lTemp(0.00);
        for ( int i = 0; i != 32; ++i ) {
          lTemp += mSizeEta[i];
          mMappingEta[i] = lTemp; //mMappingEta is the sum over all the divisions up to this point
          //std::cout << " " << mSizeEta[i] << "\t" << mMappingEta[i] << std::endl; //debug info to test we end up at eta=5
        }

      }

    double eta( const int& iEta ) {
      int lIndex = abs( iEta )-1;
      double eta = mMappingEta[lIndex] - ( 0.5*mSizeEta[lIndex] );

      if ( iEta > 0 ) {
        return eta;
      } else {
        return -eta;
      }
    }

    double phi( const int& iPhi ) {
      //return ( double( iPhi )-0.5 ) * mPhiSize;

      //return bin-centre
      double phi_degrees=0.0;
      if(iPhi < 37) { //i.e. 1-36
        phi_degrees = (5.0*(double)iPhi) - 2.5;
      } else { //i.e. 37 - 72
        phi_degrees = -1.0 * (((73.0 - (double)iPhi) * 5.0) - 2.5);
      }

      return (phi_degrees / 180.0) * TMath::Pi(); 

    }

    int iEta( const double & eta) {

      for(unsigned int i=0; i<32; i++) {
        //std::cout << fabs(eta) << ", " << mMappingEta[i] << std::endl;
        if(mMappingEta[i] > fabs(eta)) {
          if(eta < 0) { return -1 * (i+1); }
          else { return (i+1); }
        }
      }
      return 0;
    }

    int iPhi( const double & phi) {
      //assumes phi is in radians and is -pi to +pi
      //first, convert to degrees
      double phi_degrees = phi *  180.0 / TMath::Pi();
      int iPhi=0;
      if(phi_degrees < 0) {
        iPhi = 72 - (int)(fabs(phi_degrees) / 5.0);
      } else {
        iPhi = 1 + (int)(phi_degrees / 5.0);
      } 
      return iPhi;
    }

    int new_iEta(const int & iEta) { //converts from iEta to new_iEta (-28,0,+28 --> 0-56)
      int newEtaIndex1 = iEta + 27;
      if ( iEta < 0 ) { newEtaIndex1++; }
      return newEtaIndex1;
    }

    int new_iPhi(const int & iPhi) { //converts from iPhi to new_iPhi
      return iPhi - 1;
    }

    int old_iEta(const int & new_iEta) { //converts from new_iEta to iEta (0-56 --> -28,0,+28)
      int oldiEta = new_iEta - 27;
      if (oldiEta <= 0) { oldiEta--; }
      return oldiEta;
    }

    int old_iPhi(const int & new_iPhi) { //converts from new_iPhi to iPhi
      return new_iPhi + 1;
    }

    double towerEtaSize( const int& iEta ) {
      return mSizeEta[ abs( iEta )-1 ];
    }

    double towerPhiSize( const int& iPhi ) {
      return mPhiSize;
    }

    int gct_iPhi(const int& iPhi){ //Returns the gct index for phi (1-72 --> 0-17)
      
      if(iPhi > 70 || iPhi < 3) return 0;
      else return int(double(iPhi+1.0)/4.0);

    }

    int gct_iEta(const int& iEta){ //Returns the gct index for eta (-28,0,28 --> 4-17)

      int niEta = new_iEta(iEta);
      return int(double(niEta)/4.0) + 4;

    }
};

#endif

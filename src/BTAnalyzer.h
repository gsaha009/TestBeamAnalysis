#ifndef BTAnalyzer_h
#define BTAnalyzer_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1.h>
#include <TFile.h>
#include <TString.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

// Header file for the classes stored in the TTree if any.
#include<iostream>
#include<map>
#include "Cluster.h"
#include "Stub.h"
#include "Track.h"
#include "ModuleHistos.h"

//Even == Top
//Odd == Bottom
enum SENSORID{
  Mod1TOP    =30, 
  Mod1BOTTOM =31, 
  Mod2TOP    =32, 
  Mod2BOTTOM =33
};

class BTAnalyzer {
public :
  BTAnalyzer(TString inFile, unsigned int nMods, TString sensorType, TString outFile);
   virtual ~BTAnalyzer();
   bool Init();
   void Loop();
   void offlineclusterizer(std::vector<int>& hits, const unsigned int nCbc, const unsigned int nStripsPerCBC, std::vector<cluster>& clusVec );
   void stubSimulator(const std::vector<cluster>& seeding, const std::vector<cluster>& matching, 
                      std::vector<stub>& stubVec, const unsigned int clswCut = 4, const float window = 5.);
   void SaveHistos();
   void Reset();

   
private :
   void bookHistogram();
   void drawHistogram();   
   
   TH1I* topSensorHitProfile;
   TH1I* bottomSensorHitProfile;   
   
   TString inFile_;
   TString outFile_;
   TTreeReader* treeReader_;     //!pointer to the analyzed DUT TTree or TChain
   TTreeReader* tktreeReader_;   //!pointer to the analyzed DUT TTree or TChain
   TTreeReader* fptreeReader_;   //!pointer to the analyzed DUT TTree or TChain
   TTreeReader* cbcstubtreeReader_;   //!pointer to the analyzed DUT TTree or TChain

   unsigned int nmods_;
   TString sensorType_;

   TTreeReaderValue< std::vector< int > >* row_;//which side of the module?
   TTreeReaderValue< std::vector< int > >* col_;//which channel
   TTreeReaderValue< std::vector< int > >* iden_;//identifies top and bottom sensor#top=30#bottom=31
   
   //reference plane
   TTreeReaderValue< std::vector< double > >  *xPos_fp;
   TTreeReaderValue< std::vector< double > >  *yPos_fp;
   TTreeReaderValue< std::vector< int > >     *sensorId_fp;
   TTreeReaderValue< std::vector< int > >    *clustersize_fp;


   ModuleHistos*  mhists_[4];
   TrackHistos*   tkhists_;

   TTreeReaderValue< std::vector< double > >* xPos_;
   TTreeReaderValue< std::vector< double > >* yPos_;
   TTreeReaderValue< std::vector< double > >* dxdz_;
   TTreeReaderValue< std::vector< double > >* dydz_;
   TTreeReaderValue< std::vector< double > >* chi2_;
   TTreeReaderValue< std::vector< double > >* ndof_;
   TTreeReaderValue< std::vector< int > >* tkiden_;
   TTreeReaderValue< std::vector< int > >* trackNum_;

   TTreeReaderValue< std::vector< int > > *cbcstubPos_;
   TTreeReaderValue< std::vector< int > > *Fe_;
   TTreeReaderValue< std::vector< int > > *Cbc_;
   TTreeReaderValue< std::vector< int > > *stubBend_;

   bool raneventloop_;
};

#endif


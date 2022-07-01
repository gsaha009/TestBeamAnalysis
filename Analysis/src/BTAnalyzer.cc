#include "BTAnalyzer.h"
#include <TH2.h>
#include <TH1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <sstream>
#include <cstdlib>
#include "TString.h"

//#define DEBUG
//#define DEBUG2
//#define DEBUG3
//#define DEBUG4

//BTAnalyzer::BTAnalyzer(TString inFile, TString sensorType, TString outFile) : treeReader_(0) {
BTAnalyzer::BTAnalyzer(TString inFile, TString outFile) : treeReader_(0) {
  inFile_ = inFile;
  //sensorType_ = sensorType;
  outFile_ = outFile;
  raneventloop_ = false;
  if(!Init()) {
    std::cout << "Initialization error\n";
    std::exit(EXIT_FAILURE);
  } 
}



BTAnalyzer::~BTAnalyzer() {
  std::cout << "Destroying pointers!!" << std::endl;
  if (treeReader_) delete treeReader_;
  if (tktreeReader_) delete tktreeReader_;
  if (fptreeReader_) delete fptreeReader_;
  if (cbcstubtreeReader_) delete cbcstubtreeReader_;
  //if (mhists_) delete mhists_;
  //if (tkhists_) delete tkhists_;
}




// ---------------------- Initialize all the necessary variables ----------------------- //
bool BTAnalyzer::Init() {
  TFile* fin = TFile::Open(inFile_);
  if(!fin)   {
    std::cout << "File doesn't exist\n";
    std::exit(EXIT_FAILURE);
  }

  //read DUT info
  treeReader_ = new TTreeReader("rawdata",fin);
  if(!treeReader_) {
    std::cout << "DUT tree not found in file\n";
    std::exit(EXIT_FAILURE);
  }

  //tkhists
  tkhists_ = new TrackHistos();

  //loop over modules
  std::cout << "Reading branch for module: 0" << std::endl;
  row_  = new TTreeReaderValue< std::vector< int > >(*treeReader_, "row");
  col_  = new TTreeReaderValue< std::vector< int > >(*treeReader_, "col");
  iden_ = new TTreeReaderValue< std::vector< int > >(*treeReader_, "iden");

  //create the container of histograms for this module
  std::string temp = "Module_" + std::to_string(0);
  mhists_ = new ModuleHistos(temp.c_str());
  mhists_->bookHistos();

  //read tracking info
  tktreeReader_ = new TTreeReader("tracks",fin);
  if(!tktreeReader_) {
    std::cout << "Track tree not found in file\n";
    std::exit(EXIT_FAILURE);
  }
  xPos_     = new TTreeReaderValue< std::vector< double > >(*tktreeReader_, "xPos");
  yPos_     = new TTreeReaderValue< std::vector< double > >(*tktreeReader_, "yPos");
  dxdz_     = new TTreeReaderValue< std::vector< double > >(*tktreeReader_, "dxdz");
  dydz_     = new TTreeReaderValue< std::vector< double > >(*tktreeReader_, "dydz");
  chi2_     = new TTreeReaderValue< std::vector< double > >(*tktreeReader_, "chi2");
  ndof_     = new TTreeReaderValue< std::vector< double > >(*tktreeReader_, "ndof");
  tkiden_   = new TTreeReaderValue< std::vector< int > >(*tktreeReader_, "iden");
  trackNum_ = new TTreeReaderValue< std::vector< int > >(*tktreeReader_, "trackNum");

  // read tracking fit points info
  fptreeReader_ = new TTreeReader("fitpoints",fin);   
  if(!fptreeReader_) {
    std::cout << "Track tree not found in file\n";
    std::exit(EXIT_FAILURE);
  }
  xPos_fp             = new TTreeReaderValue< std::vector< double > >(*fptreeReader_, "xPos");
  yPos_fp             = new TTreeReaderValue< std::vector< double > >(*fptreeReader_, "yPos");
  sensorId_fp         = new TTreeReaderValue< std::vector< int > >(*fptreeReader_, "sensorId");
  clustersize_fp      = new TTreeReaderValue< std::vector< int > >(*fptreeReader_, "clustersize");

  //read cbc stubinfo
  cbcstubtreeReader_  = new TTreeReader("stubs",fin);
  if(!cbcstubtreeReader_) {
    std::cout << "cbc stub info not found in file\n";
  }
  //ncbcStubs_          = new TTreeReaderValue<Int_t>(*cbcstubtreeReader_, "nStubs");
  cbcstubPos_         = new TTreeReaderValue< std::vector< int > >(*cbcstubtreeReader_, "stubPos");
  Fe_                 = new TTreeReaderValue< std::vector< int > >(*cbcstubtreeReader_, "Fe");
  Cbc_                = new TTreeReaderValue< std::vector< int > >(*cbcstubtreeReader_, "Cbc");
  stubBend_           = new TTreeReaderValue< std::vector< int > >(*cbcstubtreeReader_, "stubBend");
  std::cout<<"Initialized...\n";
  return true;
}




//clear all the containers 
void BTAnalyzer::Reset() {
  //loop over modules
}



// --------------------------------------- event loop function ---------------------------------------------- //
void BTAnalyzer::Loop()
{
  TTreeReaderValue<Int_t> ncbcStubs_(*cbcstubtreeReader_, "nStubs");
  TTreeReaderValue<Int_t> tdc_(*treeReader_, "TDC");
  
  unsigned int noTk               = 0;
  unsigned int tracksPerEvent     = 0;
  unsigned int reftracksPerEvent  = 0;
  unsigned int nRefTracks         = 0;
  unsigned int nRefTracksTop         = 0;
  unsigned int nRefTracksBot         = 0;
  unsigned int nMatchedTracks     = 0;
  unsigned int nMatchedTracksTop     = 0;
  unsigned int nMatchedTracksBot     = 0;
  unsigned int nMatchedcbcStubs   = 0;
  unsigned int nMatchedcbcStubsTop   = 0;
  unsigned int nMatchedcbcStubsBot   = 0;
  unsigned int nMatchedrecoStubs  = 0;
  unsigned int nMatchedrecoStubsTop  = 0;
  unsigned int nMatchedrecoStubsBot  = 0;

  //test
  std::vector<float> clsBotFe0_allevt;
  std::vector<float> clsBotFe1_allevt;
  int evidx = -1;
  // Main Event Loop
  while (treeReader_->Next() && tktreeReader_->Next() && fptreeReader_->Next() && cbcstubtreeReader_->Next()) {
    evidx++;
    Reset();
    int tdc = *(tdc_.Get());
    std::vector<Track>  cleanedTk;
    std::vector<int> ref_tracks;
    tkhists_->nTk->Fill((xPos_->Get())->size());
#ifdef DEBUG
    std::cout<<"\n\n\nEvents : "<<evidx<<"\t"
	     <<"No of Tracks : "<<(xPos_->Get())->size()<<"\n";
#endif
    noTk += (xPos_->Get())->size(); // to get the total number of tracks saved in tree

    // Track information. Match tracks at fei4 first. Then use only selected tracks for analysis
    for (unsigned int it = 0; it < (tkiden_->Get())->size() ; it++)  {
      int identk = (tkiden_->Get())->at(it);
#ifdef DEBUG
      std::cout<<"  Track no : "<<it<<"\t"<<"identity of track : "<<identk<<"\n";
#endif
      if(identk != 20)   continue; // 8, 20
      double xtk = (xPos_->Get())->at(it);
      double ytk = (yPos_->Get())->at(it);
      tracksPerEvent++;
    
      // Loop over ref hits i.e. sensor ids of the fit points
#ifdef DEBUG
      std::cout<<"  Loop over ref hits ...\n";
#endif
      for (unsigned int iref = 0; iref < (sensorId_fp->Get())->size() ; iref++)  {
        int refId = (sensorId_fp->Get())->at(iref);
#ifdef DEBUG
	std::cout<<"    SensorID of corresponding fit point that should match the sensorId of FeI4 : "<<refId<<std::endl;
#endif
	if(refId != 20)    continue; //sensorId 10 is FeI4
        double xref = (xPos_fp->Get())->at(iref);
        double yref = (yPos_fp->Get())->at(iref);
        //residual fei4
        double xres = xref - xtk;
        double yres = yref - ytk;
#ifdef DEBUG
	std::cout<<"    Residuals : \n";
	std::cout<<"    X_res : "<<xres<<"\t"<<"Y_res : "<<yres<<std::endl;
#endif
	tkhists_->rawResXfei4->Fill(xres);
	tkhists_->rawResYfei4->Fill(yres);
        if (std::fabs(xres) < 1. && std::fabs(yres) < 1.) {
          tkhists_->resXfei4->Fill(xres);
          tkhists_->resYfei4->Fill(yres);
          tkhists_->refx_x->Fill(xtk, xres);
          tkhists_->refx_y->Fill(xtk, yres);
          tkhists_->refy_x->Fill(ytk, xres);
          tkhists_->refy_y->Fill(ytk, yres);
        }
        //residual cut at FeI4  
        if(std::fabs(xres) > 0.25 or std::fabs(yres) > 0.1 )   continue;
#ifdef DEBUG
	std::cout<<"    Trk FOUND !\n";
#endif
	//if(std::fabs(xres) > 0.20 or std::fabs(yres) > 0.08 )   continue;
        tkhists_->matchedRefHitMap->Fill(xtk, ytk);
        ref_tracks.push_back((trackNum_->Get())->at(it));
        reftracksPerEvent++;
	cleanedTk.push_back(Track((xPos_->Get())->at(it),
				  (yPos_->Get())->at(it),
				  (dxdz_->Get())->at(it),
				  (dydz_->Get())->at(it),
				  (chi2_->Get())->at(it),
				  (ndof_->Get())->at(it),
				  (tkiden_->Get())->at(it)));
      }
    } 
    // Loop over tracks ends here and 
    // we found the number of tracks that  
    tkhists_->ncleanedTk->Fill(cleanedTk.size());

    ModuleHistos& mh = *(mhists_);
    //create temp containers of hits per sensor and fe
    std::vector<int> hitsBotfe0;
    std::vector<int> hitsBotfe1;
    std::vector<int> hitsTopfe0;
    std::vector<int> hitsTopfe1;
    std::vector<cluster> clsBotfe0;
    std::vector<cluster> clsBotfe1;
    std::vector<cluster> clsTopfe0;
    std::vector<cluster> clsTopfe1;
    // test
    std::vector<cluster> clsBot;
    std::vector<cluster> clsTop;
    
    std::vector<stub> cbcstubfe0;
    std::vector<stub> cbcstubfe1;
    std::vector<stub> recostubfe0;
    std::vector<stub> recostubfe1;
    
    // TOP=30, BOTTOM=31 [Module: 1]
    // row 1 :: FE0
    // row 0 :: FE1
    // for fe0, cbc strip = 1016 - strip
#ifdef DEBUG2
    (ref_tracks.size() > 0 ) ? std::cout<<"  "<<ref_tracks.size()<<" RefTracks FOUND ! \n\n" : std::cout<<"  NO RefTrack FOUND ! \n\n"; 
    std::cout<<"  Ready to prepare offline clusters !!!\n" ;
    std::cout<<"  nColumns (/ nChannels fired): "<<(col_->Get())->size()<<"\n";
    if ((col_->Get())->size() > 0) std::cout<<"  Loop over hits \n";
#endif
    for (unsigned int nh = 0; nh < (col_->Get())->size() ; nh++)  {
      const int sid =  (iden_->Get())->at(nh);
      const int r = (row_->Get())->at(nh);  // Row is either 0 or 1
      const int ch = (col_->Get())->at(nh); // Channel is basically the cbc channel has signal. Column means all the channels having hit
      
#ifdef DEBUG2
      std::cout<<"    SensorId [30:Top, 31:Bot]: "<<sid
	       <<"    Row [0:Fe1, 1:Fe0]: "<<r
	       <<"    Channel : "<<ch
	       <<std::endl;
#endif
      //TOP=30, BOTTOM=31
      //For bottom sensor, mask channel 794 for row==0 [for downstream module only??]
      //if( sid == SENSORID::Mod1TOP) {
      if (sid == 30) {
	if(r == 1) hitsTopfe0.push_back(ch);
	else if (r == 0) hitsTopfe1.push_back(ch);
      }
      //else if( sid == SENSORID::Mod1BOTTOM) { 
      if( sid == 31) { 
	if(r == 1) hitsBotfe0.push_back(ch);
	else if (r == 0 and ch != 794) hitsBotfe1.push_back(ch);
      }
    }//finish loop over hits
    
#ifdef DEBUG2
    std::cout<<"  nHitsBotFe0 : "<<hitsBotfe0.size()<<"\t"
	     <<"  nHitsBotFe1 : "<<hitsBotfe1.size()<<"\n"
	     <<"  nHitsTopFe0 : "<<hitsTopfe0.size()<<"\t"
	     <<"  nHitsTopFe1 : "<<hitsTopfe1.size()
	     <<std::endl;
#endif
    //Plot hit prop
    mh.bottomS.nhits_Fe0->Fill(hitsBotfe0.size());
    mh.bottomS.nhits_Fe1->Fill(hitsBotfe1.size());
    mh.topS.nhits_Fe0->Fill(hitsTopfe0.size());
    mh.topS.nhits_Fe1->Fill(hitsTopfe1.size());
    for (auto& h: hitsBotfe0) mh.bottomS.hitpos_Fe0->Fill(h);
    for (auto& h: hitsBotfe1) mh.bottomS.hitpos_Fe1->Fill(h);
    for (auto& h: hitsTopfe0) mh.topS.hitpos_Fe0->Fill(h);
    for (auto& h: hitsTopfe1) mh.topS.hitpos_Fe1->Fill(h);
    
    //create clusters
    offlineclusterizer(hitsBotfe0, 8, 127, clsBotfe0);
    offlineclusterizer(hitsBotfe1, 8, 127, clsBotfe1);
    offlineclusterizer(hitsTopfe0, 8, 127, clsTopfe0);
    offlineclusterizer(hitsTopfe1, 8, 127, clsTopfe1);

#ifdef DEBUG2
    if (hitsBotfe0.size() + hitsBotfe1.size() + hitsTopfe0.size() + hitsTopfe1.size() > 0) std::cout<<"    Clusterisation done !\n";
    if (clsBotfe0.size() + clsBotfe1.size() + clsTopfe0.size() + clsTopfe1.size() > 0) {
      std::cout<<"  nClsBotFe0 : "<<clsBotfe0.size()<<"\t"
	       <<"  nClsBotFe1 : "<<clsBotfe1.size()<<"\n"
	       <<"  nClsTopFe0 : "<<clsTopfe0.size()<<"\t"
	       <<"  nClsTopFe1 : "<<clsTopfe1.size()
	       <<std::endl;
    }
#endif

    //test
    for (auto& cls : clsBotfe0) {
      clsBot.push_back(cls);
      mh.bottomS.clsposfe0   ->Fill(cls.center());
      mh.bottomS.clspos      ->Fill(cls.center());
      mh.bottomS.clswidthfe0 ->Fill(cls.size());
#ifdef DEBUG2
      std::cout<<"    fstrip : "<< cls.firstStrip() <<"\t"
	       <<"    center : "<< cls.center()     <<"\t"
	       <<"    size   : "<< cls.size()       <<"\t"
	       <<"    column : "<< cls.column()
	       <<std::endl;
#endif
    }
    for (auto& cls : clsBotfe1) {
      int fstrip = cls.firstStrip() + 1016;
      cluster clust(fstrip, cls.column(), cls.size());
      clsBot.push_back(clust);
      mh.bottomS.clsposfe1   ->Fill(cls.center());
      mh.bottomS.clspos      ->Fill(cls.center());
      mh.bottomS.clswidthfe1 ->Fill(cls.size());
#ifdef DEBUG2
      std::cout<<"    fstrip : "<< clust.firstStrip() <<"\t"
	       <<"    center : "<< clust.center()     <<"\t"
	       <<"    size   : "<< clust.size()       <<"\t"
	       <<"    column : "<< clust.column()
	       <<std::endl;
#endif
    }
    for (auto& cls : clsTopfe0) {
      clsTop.push_back(cls);
      mh.topS.clsposfe0     ->Fill(cls.center());
      mh.topS.clspos        ->Fill(cls.center());
      mh.topS.clswidthfe0   ->Fill(cls.size());
#ifdef DEBUG2
      std::cout<<"    fstrip : "<< cls.firstStrip() <<"\t"
	       <<"    center : "<< cls.center()     <<"\t"
	       <<"    size   : "<< cls.size()       <<"\t"
	       <<"    column : "<< cls.column()
	       <<std::endl;
#endif
    }
    for (auto& cls : clsTopfe1) {
      int fstrip = cls.firstStrip() + 1016;
      cluster clust(fstrip, cls.column(), cls.size());
      clsTop.push_back(clust);
      mh.topS.clsposfe1     ->Fill(cls.center());
      mh.topS.clspos        ->Fill(cls.center());
      mh.topS.clswidthfe1   ->Fill(cls.size());
#ifdef DEBUG2
      std::cout<<"    fstrip : "<< clust.firstStrip() <<"\t"
	       <<"    center : "<< clust.center()     <<"\t"
	       <<"    size   : "<< clust.size()       <<"\t"
	       <<"    column : "<< clust.column()
	       <<std::endl;
#endif
    }
    // correlation
    for(auto& h0: clsBotfe0) {
      for(auto& h1: clsTopfe0)
	mh.clusterPoscorrfe0->Fill(h0.center(), h1.center());
    }
    for(auto& h0: clsBotfe1) {
      for(auto& h1: clsTopfe1)
	mh.clusterPoscorrfe1->Fill(h0.center(), h1.center());
    }
    for(auto& h0: clsBot) {
      for(auto& h1: clsTop)
	mh.clusterPoscorr->Fill(h0.center(), h1.center());
    }

#ifdef DEBUG3
    std::cout<<"\n  Ready to prepare CBC stubs !\n";
    std::cout<<"  nCBC_Stubs : "<<(cbcstubPos_->Get())->size()<<"\n";
#endif
    //CBC Stub information
    for(unsigned int is = 0; is <  (cbcstubPos_->Get())->size(); is++) {
      int fe   = (Fe_->Get())->at(is);
      int cbc  = (Cbc_->Get())->at(is);
      int pos  = (cbcstubPos_->Get())->at(is);
      int bend = (stubBend_->Get())->at(is);
      float strip = cbc*127 + pos/2.;
      if (fe == 0) {
	strip = 1016 - strip;
	cbcstubfe0.emplace_back(stub(strip, fe, bend));
      } 
      else if (fe == 1) 
	cbcstubfe1.emplace_back(stub(strip, fe, bend));
    }//cbc stub loop
    

    //plot cbc stub prop
    mh.ncbcstubsfe0->Fill(cbcstubfe0.size());
    std::vector<float>cbcStubPos;
    std::vector<float>cbcStubBend;
    for(auto& s : cbcstubfe0) {
      //float stubPos = (s.center()-508.5)*0.09;
      cbcStubPos.push_back(/*stubPos*/s.realPos());
      cbcStubBend.push_back(s.bend());
      mh.cbcStubPos -> Fill (/*stubPos*/s.realPos());
      mh.cbcStubposfe0->Fill(s.center());
      //std::cout<<"fe0 : "<<s.bend()<<"\n";
      mh.cbcStubbendfe0->Fill(s.bend());
#ifdef DEBUG3
      std::cout<<"    CBC Stubs info at Fe0 "<<"\n"
	       <<"      posFe0   : "<<s.center()  <<"\t"
	       <<"      realPos  : "<<s.realPos() <<"\t"
	       <<"      bend     : "<<s.bend()
	       <<std::endl;
#endif
      for (auto& h: clsBotfe0) {
	mh.clsvsStubPoscorrfe0->Fill(h.center(), s.center());
      }
    }
    mh.ncbcstubsfe1->Fill(cbcstubfe1.size());
    for(auto& s : cbcstubfe1) {
      //float stubPos = (s.center()-508.5)*0.09;
      cbcStubPos.push_back(/*stubPos*/s.realPos());
      cbcStubBend.push_back(s.bend());
      mh.cbcStubPos -> Fill (/*stubPos*/s.realPos());
      mh.cbcStubposfe1->Fill(s.center());
      //std::cout<<"fe1 : "<<s.bend()<<"\n";
      mh.cbcStubbendfe1->Fill(s.bend());
#ifdef DEBUG3
      std::cout<<"    CBC Stubs info at Fe1 "<<"\n"
	       <<"      posFe1   : "<<s.center()  <<"\t"
	       <<"      realPos  : "<<s.realPos() <<"\t"
	       <<"      bend     : "<<s.bend()
	       <<std::endl;
#endif
      for (auto& h: clsBotfe1) {
	mh.clsvsStubPoscorrfe1->Fill(h.center(), s.center());
      }
    }
    
    //now do offline stub reco 
    //seeding, matching, stubVec

#ifdef DEBUG3
    std::cout<<"  Offline stub simulation\n";
#endif
    stubSimulator(clsBotfe0, clsTopfe0, recostubfe0);//bottom is seeding
    stubSimulator(clsBotfe1, clsTopfe1, recostubfe1);//bottom is seeding
    
    //plot reco stub prop
    mh.nrecostubsfe0->Fill(recostubfe0.size());
    std::vector<float>recoStubPos;
    for(auto& rs : recostubfe0) {
      //float stubPos = (rs.center()-508.5)*0.09;
      recoStubPos.push_back(rs.realPos());
      mh.recoStubPos -> Fill (rs.realPos());
      mh.recoStubposfe0->Fill(rs.center());
#ifdef DEBUG3
      std::cout<<"    Stub info Fe0 .. \n"
	       <<"      col     : "<<rs.column()
	       <<"      bend    : "<<rs.bend()
	       <<"      center  : "<<rs.center()
	       <<"      realpos : "<<rs.realPos()
	       <<std::endl;
#endif
      for (auto& cs: cbcstubfe0) {
	mh.rsvscsPoscorrfe0->Fill(cs.center(), rs.center());
      }
    }
    mh.nrecostubsfe1->Fill(recostubfe1.size());
    
    for(auto& rs : recostubfe1) {
      float stubPos= 0.09 + (rs.center()-508.5)*0.09; 
      mh.recoStubPos -> Fill (stubPos);
      recoStubPos.push_back(stubPos);
      mh.recoStubposfe1->Fill(rs.center());
#ifdef DEBUG3
      std::cout<<"    Stub info Fe1 .. \n"
	       <<"      col     : "<<rs.column()
	       <<"      bend    : "<<rs.bend()
	       <<"      center  : "<<rs.center()
	       <<"      realpos : "<<rs.realPos()
	       <<std::endl;
#endif
      for (auto& cs: cbcstubfe1) {
	mh.rsvscsPoscorrfe1->Fill( cs.center(), rs.center() );
      }
    }
    
#ifdef DEBUG3
    std::cout<<"\n  TDC >>----> "<<tdc<<std::endl;
#endif
    // now Tk to DUT matching
    // If there is ref hit, take the last one
    // if (ref_tracks>-1)
    bool tdcok {false};
    //tdcok = true;
    //tdcok = tdc == 6; // new beamtest
    //tdcok = tdc == 2 or tdc == 3 or tdc == 4; // old beamtest
    tdcok = true;
    if(!tdcok)   continue;

    /*
    int senID = 0;
    //std::cout<<sensorType_<<"\t"<<senID<<"\n";
    if (sensorType_ == "bottom")   senID = 31;
    else if (sensorType_ == "top") senID = 30;
    else {
      std::cout<<"sensortype has not been mentioned correctly! ['bottom' || 'top']\n";
      senID = 31; //bottom by default
    }
    */

    std::vector<int>senIds {31,30};
#ifdef DEBUG4
    std::cout<<"  Now iterate over the two sensors 30 [top] and 31 [bottom] ..."<<std::endl;
#endif
    for (int senId : senIds) {
      if (senId == 31) sensorType_ = "bottom";
      else if (senId == 30) sensorType_ = "top";
#ifdef DEBUG4
      std::cout<<"    senId   : "<<senId<<"\t"
	       <<"    senType : "<<sensorType_
	       <<std::endl;
#endif
      int nCls = 0;
      bool istop {false};
      bool isbot {false};
      std::vector<double>refTkPos;
    
      // Loop over reference tracks
#ifdef DEBUG4
      std::cout<<"    Loop over ref tracks ... \n";
#endif
      for (unsigned int irt =0; irt < ref_tracks.size(); irt++) {
	int ref_track = ref_tracks[irt]; // get ref track index
	int trackID = -1;
#ifdef DEBUG4
	std::cout<<"      Ref trk no : "<<ref_track<<"\n"
		 <<"      Checking if it has hit on DUT or not, so loop over total tracks to find the same ref trk and check its id\n";
#endif      
	// Loop over total tracks to check if a track has DUT hit or not
	for (unsigned int it=0; it < (tkiden_->Get())->size(); it++) {
	  // Find the ref track index from the track collection and check if that track 
	  // has a hit at DUT lower or upper sensor
#ifdef DEBUG4
	std::cout<<"        trk no   : "<<(trackNum_->Get())->at(it)<<"\t"
		 <<"        trk iden : "<<(tkiden_->Get())->at(it)
		 <<std::endl;
#endif
	  if ((trackNum_->Get())->at(it) == ref_track && (tkiden_->Get())->at(it) == senId) {
	    trackID=it;
	    nRefTracks++;
	    if (senId == 30) nRefTracksTop++;
	    else if (senId == 31) nRefTracksBot++;
#ifdef DEBUG4
	    std::cout<<"        ..matched : "<<"trk no : "<<(trackNum_->Get())->at(trackID)<<"\t"<<"senId : "<<senId<<"\n";
#endif
	    break;
	  }
	}
	//Track position on DUT plane
	double txpos = (xPos_->Get())->at(trackID);
	double typos = (yPos_->Get())->at(trackID);
#ifdef DEBUG4
	std::cout<<"      trk x-pos on DUT : "<<txpos<<"\t"
		 <<"      trk y-pos on DUT : "<<typos
		 <<std::endl;
#endif
	mh.tkatDUTmap->Fill(txpos, typos);
	mh.tkxatDUT->Fill(txpos);
	mh.tkyatDUT->Fill(typos);
	//Fill denominator tdc histogram
	mh.htdc->Fill(tdc);
	bool isMatched = false;
	bool tkMcount[10] = {0};
	
	// now loop over sensor hits to get the matched hit with the track selected
	// The sensor hits have info of DUT sensors and FeI4
	// So it is necessary to ensure that the sensorID must be of DUT
#ifdef DEBUG4
	std::cout<<"      Now loop over sensor hits to get the matched hit with the track selected"<<"\n"
		 <<"      The sensor hits have info of DUT sensors and FeI4"<<"\n"
		 <<"      So it is necessary to ensure that the sensorID must be of DUT\n";
#endif
	for (unsigned int iref = 0; iref < (sensorId_fp->Get())->size() ; iref++)  {
	  int refId = (sensorId_fp->Get())->at(iref);
	  if(refId != senId)    continue; //sensorId :: 30,31 for mod1 || 32,33 for mod2 :: is DUT
	  // Get the track fit points on DUT
	  double xdut = (xPos_fp->Get())->at(iref);
	  double ydut = (yPos_fp->Get())->at(iref);
	  unsigned int clsize  = (clustersize_fp->Get())->at(iref);
	  //residual DUT
	  double xres = xdut - txpos;
	  double yres = ydut - typos;
#ifdef DEBUG4
	  std::cout<<"        trk fit point x-pos on DUT : "<<xdut<<"\t"
		   <<"        trk fit point y-pos on DUT : "<<ydut<<"\n"
		   <<"        Now choose nearest strip approach to get the closest fit-point of trk ... \n";

#endif
	  //double minRes = xres; // extra line to keep rest of the code the same 
	  // --------------------------------------- NSA ----------------------------------------- //
	  double ce = xdut - 0.5*clsize*0.09; //Reached at the left edge of the cluster
	  double minRes = 999.9;
	  //Finding the strip positions for this cluster
	  for (size_t is = 0; is < clsize; ++is){
	    double pos = ce + (is + 0.5) * 0.09; 
	    double res = std::fabs(pos - txpos);
	    if (res < minRes) minRes = res;
	  }
#ifdef DEBUG4
	  std::cout<<"          Minimum residual : "<<minRes<<std::endl;
#endif
	  // ------------------------------------------------------------------------------------- //
	  mh.tkVsDUTxcorr->Fill(txpos, xdut);
	  mh.tkVsDUTycorr->Fill(typos, ydut);
	  mh.resXdut->Fill(xres);
	  mh.resYdut->Fill(yres);
	  for(unsigned int ires = 0; ires < 10; ires++) {
	    float cut = (ires+1) * 0.05;
	    if ( minRes < cut) 
	      tkMcount[ires] = true;
	  }
	  //This cut has to be varied to plot efficiency as a function of cut
	  //if ( minRes < 0.15 /*&& fabs(dy)<2.5*/) {
	  if ( minRes < 0.20 /*&& fabs(dy)<2.5*/) {
	    nMatchedTracks++;
	    mh.nMatchedTk -> Fill(0);
	    nCls++;
	    isMatched = true;
	    if (sensorType_ == "top") {
	      nMatchedTracksTop++;
	      mh.topS.clswidth_alltdc->Fill(clsize);
	      if (tdc == 2) mh.topS.clswidth_tdc2->Fill(clsize);
	      else if (tdc == 3) mh.topS.clswidth_tdc3->Fill(clsize);
	      else if (tdc == 4) mh.topS.clswidth_tdc4->Fill(clsize);
	      istop = true;
	    }
	    else if (sensorType_ == "bottom") {
	      nMatchedTracksBot++;
	      mh.bottomS.clswidth_alltdc->Fill(clsize);
	      if (tdc == 2) mh.bottomS.clswidth_tdc2->Fill(clsize);
	      else if (tdc == 3) mh.bottomS.clswidth_tdc3->Fill(clsize);
	      else if (tdc == 4) mh.bottomS.clswidth_tdc4->Fill(clsize);
	      isbot = true;
	    }
	  }
	}//end loop over sensor hits
	for(unsigned int ires = 0; ires < 10; ires++) 
	  if(tkMcount[ires]) mh.tkmatchCounter->Fill((ires+1)*50);
	
	if(isMatched)   {
	  mh.htdc_matched->Fill(tdc);
	  refTkPos.push_back(txpos);
	  mh.tkxatDUT_matched->Fill(txpos);
	  mh.tkyatDUT_matched->Fill(typos);
	  
	  
	  // matched cbc Stubs
	  double min = 999.9;
	  double cbcStPos = 0.0;
	  float cbcStBend = 19.0;
	  bool hasMatCbc = false;
	  bool hasMatRec = false;
	  
	  int index = 0;
	  for (auto& s: cbcStubPos) {
	    index++;
	    double diff = s - txpos;
	    if (std::abs(diff) < std::abs(min)) {
	      cbcStPos = s;
	      cbcStBend = cbcStubBend[index];
	      min = diff;
	    }
	  }
	  mh.cbcStubTrackPosDiff->Fill(min);
	  if (std::abs(min) < 0.20/* && std::abs(cbcStBend) < 2*/) {
	    mh.tkcbcStubRes->Fill(min);
	    mh.matchedCbcStubPos->Fill(cbcStPos);
	    //std::cout<<"cbcStPos: "<<cbcStPos<<"\t";
	    mh.matchedCbcStubBend->Fill(cbcStBend);
	    nMatchedcbcStubs++;
	    if (sensorType_ == "bottom") nMatchedcbcStubsBot++;
	    else if (sensorType_ == "top") nMatchedcbcStubsTop++;
	    mh.nMatchedcbcSt->Fill(0);
	    hasMatCbc = true;
	  }
	  // matched reco Stubs
	  min = 999.9;
	  double recoStPos = 0.0;
	  for (auto& s: recoStubPos) {
	    double diff = s - txpos;
	    if (std::abs(diff) < std::abs(min)) {
	      recoStPos = s;
	      min = diff;
	    }
	  }
	  mh.recoStubTrackPosDiff->Fill(min);
	  if (std::abs(min) < 0.20) {
	    mh.tkrecoStubRes->Fill(min);
	    mh.matchedRecoStubPos->Fill(recoStPos);
	    //std::cout<<"recoStPos: "<<recoStPos<<"\n";
	    nMatchedrecoStubs++;
	    if (sensorType_ == "bottom") nMatchedrecoStubsBot++;
	    else if (sensorType_ == "top") nMatchedrecoStubsTop++;
	    mh.nMatchedrecoSt->Fill(0);
	    hasMatRec = true;
	  }
	  //std::cout<<"\n";
	  if (hasMatRec && hasMatCbc) mh.stubPosCorrl->Fill(cbcStPos, recoStPos);
	}//if (isMatched) scope ends
      }//loop over ref tracks ends
    
      mh.cumlEffCbc->Fill(nMatchedTracks, float(nMatchedcbcStubs-nMatchedrecoStubs));
      mh.cumlEffReco->Fill(nMatchedTracks, float(nMatchedrecoStubs)/nMatchedTracks);
    
      if (isbot) {
	mh.bottomS.ncls_alltdc->Fill(nCls);
	if (tdc == 2) mh.bottomS.ncls_tdc2->Fill(nCls);
	else if (tdc == 3) mh.bottomS.ncls_tdc3->Fill(nCls);
	else if (tdc == 4) mh.bottomS.ncls_tdc4->Fill(nCls);
      }
      
      else if (istop) {
	mh.topS.ncls_alltdc->Fill(nCls);//new
	if (tdc == 2) mh.topS.ncls_tdc2->Fill(nCls);//new
	else if (tdc == 3) mh.topS.ncls_tdc3->Fill(nCls);//new
	else if (tdc == 4) mh.topS.ncls_tdc4->Fill(nCls);//new
      }
    }
  }//iev++;
      //event loop while end
 
  raneventloop_ = true;
  std::cout << "Total tracks =" << tracksPerEvent << std::endl; 
  std::cout << "Module [1]\n";
  std::cout << "nRefTracks     = " << nRefTracksBot << std::endl;
  std::cout << "nMatchedTracks = " << nMatchedTracksBot << std::endl;
  std::cout << "Track matching Efficiency with Clusters [NSA :: 150µ] :: Bottom sensor : " << (float)nMatchedTracksBot/nRefTracksBot*100.  << "%" << std::endl;
  std::cout << "Track matching Efficiency with Clusters [NSA :: 150µ] :: Top sensor    : " << (float)nMatchedTracksTop/nRefTracksTop*100.  << "%" << std::endl;
  std::cout << "nMatchedCBCStubs: "<< nMatchedcbcStubsBot<<"\t"<< "nMatchedRECOStubs: "<< nMatchedrecoStubsBot<<"\t"    
	    << "cbcStubEfficiency: "<<(float)nMatchedcbcStubsBot/nMatchedTracksBot*100 << "%" << "\t" 
	    << "recoStubEfficiency: "<<(float)nMatchedrecoStubsBot/nMatchedTracksBot*100 << "%" <<std::endl;  
  
}//end event loop 

//Function to create clusters from hits
void BTAnalyzer::offlineclusterizer(std::vector<int>& hits, 
				    const unsigned int nCbc, 
				    const unsigned int nStripsPerCBC, 
				    std::vector<cluster>& clusVec ) {
#ifdef DEBUG2
  if (hits.size() > 0) std::cout<<"  Making offline cluster===>>>"<<"\n";
#endif
  std::sort(hits.begin(), hits.end());
  if (hits.empty())  return; 
  unsigned int fStrip = hits.at(0);//get strip of first hit
  unsigned int col    = hits.at(0)/16 < 8 ?  0 : 1;//0-1015 is col 0 and 1016-2031 is col 1
  //what can be done? add a chip variable in cluster object and fill the cluster info here
  //should be similar to what is done in DQM/ This info can then be propagated to stubs as well.
  //std::cout<<"firstStrip: "<<fStrip<<"\t"<<"col: "<<col<<"\n";
  unsigned int size=1;
  //unsigned int edge = 8*nStripsPerCBC;
  //if (hits.size() > 1) std::cout<<"Loop over hits:::\n";
#ifdef DEBUG2
  std::cout<<"    hit: "<<0<<"\t"<<"strip: "<<fStrip<<"\t"<<"col: "<<col<<"\t"<<" for Fe1 1016 will be added later\n";
#endif

  for (unsigned int i = 1; i < hits.size(); i++){
    unsigned int icol    = hits.at(i)/16 < 8 ?  0 : 1;
#ifdef DEBUG2
    std::cout<<"    hit: "<<i<<"\t"<<"strip: "<<hits.at(i)<<"\t"<<"col: "<<icol<<"\t"<<" for Fe1 1016 will be added later\n";
#endif
    //form vectors of with hits from the same column
    //if (hits.at(i) == fStrip + size && icol == col){
    if (std::fabs(hits.at(i)-fStrip) == size && icol == col){
      size++;
    }
    else{
      cluster clust(fStrip, col, size);
      //std::cout<<"Cluster:   "<<"firstStrip: "<<fStrip<<"\t"<<"col: "<<col<<"\t"<<"clsSize: "<<size<<"\n";
      clusVec.push_back(clust);
      //reset the intial parameters
      size=1;
      fStrip = hits.at(i);//get strip of first hit
      col    = icol;
    }  
  }       
  cluster clust(fStrip, col, size);
  //std::cout<<"Cluster:   "<<"firstStrip: "<<fStrip<<"\t"<<"col: "<<col<<"\t"<<"clsSize: "<<size<<"\n";
  clusVec.push_back(clust);
}

//function to create stubs from clusters
void BTAnalyzer::stubSimulator(const std::vector<cluster>& seeding, 
			       const std::vector<cluster>& matching, 
			       std::vector<stub>& stubVec, 
                               const unsigned int clswCut, 
			       const float window )
{
  //  std::cout<<"stub simulator===>>>>\n";
  for(auto& sCls : seeding) {
    if(sCls.size() > clswCut)    continue;//cut cluster size
    //std::cout<<"sClsSize: "<<sCls.size()<<"\t"<<"sClsCenter: "<<sCls.center()<<"\n";
    for(auto& mCls :matching) {
      if(mCls.size() > clswCut)  continue;//cut cluster size
      //std::cout<<"mClsSize: "<<mCls.size()<<"\t"<<"mClsCenter: "<<mCls.center()<<"\n";
      //this is a FP operation. Wondering if it is better to save all FP positions as uint32(multipying by 2)
      //So the idea would be: a cluster at pos = 100.5 will be represented as 201
      //the window is an integer, so converting uint32 can be a better option//someone should check.
      //std::cout<<"window: "<<std::abs(sCls.center() - mCls.center())<<"\n";
      if(std::abs(sCls.center() - mCls.center()) <= window) {
	//	std::cout<<"SeedCls: "<<sCls.center()<<"\t"<<"matCls: "<<mCls.center()<<"\n";
	//        stubVec.emplace_back( stub(sCls.firstStrip(), sCls.column(), sCls.center() - mCls.center()) );
        stubVec.emplace_back( stub(sCls.center(), sCls.column(), sCls.center() - mCls.center()) );
      }
    }
  }
}

//Save histograms in outpur root file. One folder per module.
//have to think how to organize inter module correlation
void BTAnalyzer::SaveHistos() {
  if(!raneventloop_) std::cout << "You are trying to save histograms without running the event loop!" << std::endl; 
  TFile* fout = TFile::Open(outFile_, "recreate");
  fout->cd();
  tkhists_->writeHistostofile(fout, "tracks");
  fout->cd();
  mhists_->writeHistostofile(fout);
  fout->Save();
  fout->Print();
  fout->Close();
}

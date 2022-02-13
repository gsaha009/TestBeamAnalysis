#include "BTAnalyzer.h"
#include <TH2.h>
#include <TH1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <sstream>
#include <cstdlib>
#include "TString.h"

BTAnalyzer::BTAnalyzer(TString inFile, unsigned int nMods, TString sensorType, TString outFile) : treeReader_(0) {
  inFile_ = inFile;
  nmods_ = nMods;
  sensorType_ = sensorType;
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
  //tkhists_
  tkhists_ = new TrackHistos();
  //loop over modules
  for(unsigned int i = 0; i < nmods_; i++ ) {
    std::cout << "Reading branch for module:" << i << std::endl;
    row_  = new TTreeReaderValue< std::vector< int > >(*treeReader_, "row");
    col_  = new TTreeReaderValue< std::vector< int > >(*treeReader_, "col");
    iden_ = new TTreeReaderValue< std::vector< int > >(*treeReader_, "iden");
    //create the container of histograms for this module
    std::string temp = "Module_" + std::to_string(i);
    //mhists_[i] = new ModuleHistos(TString("Module_" + std::to_string(i)));
    mhists_[i] = new ModuleHistos(temp.c_str());
    mhists_[i]->bookHistos();
  }
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
  //read tracking info
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
  //ncbcStubs_ = new TTreeReaderValue<Int_t>(*cbcstubtreeReader_, "nStubs");
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
  
  unsigned int noTk = 0;
  unsigned int tracksPerEvent = 0;
  unsigned int reftracksPerEvent = 0;
  unsigned int nRefTracks = 0;
  unsigned int nMatchedTracks = 0;
  unsigned int nMatchedcbcStubs   = 0;
  unsigned int nMatchedrecoStubs  = 0;
  unsigned int nRefTracksMod1 = 0;
  unsigned int nRefTracksMod2 = 0;
  unsigned int nMatchedTracksMod1 = 0;
  unsigned int nMatchedTracksMod2 = 0;
  unsigned int nMatchedcbcStubsMod1   = 0;
  unsigned int nMatchedrecoStubsMod1  = 0;
  unsigned int nMatchedcbcStubsMod2   = 0;
  unsigned int nMatchedrecoStubsMod2  = 0;

  //test
  std::vector<float> clsBotFe0_allevt;
  std::vector<float> clsBotFe1_allevt;


  //int iev = 0;
  while (treeReader_->Next() && tktreeReader_->Next() && fptreeReader_->Next() && cbcstubtreeReader_->Next()) {
    //clear all the containers
    Reset();
    int tdc = *(tdc_.Get());

    std::vector<Track>  cleanedTk;
    std::vector<int> ref_tracks;
    tkhists_->nTk->Fill((xPos_->Get())->size());
    noTk += (xPos_->Get())->size();
    //Track information
    //match tracks at fei4 first. then use only selected tracks for analysis
    for (unsigned int it = 0; it < (tkiden_->Get())->size() ; it++)  {
      int identk = (tkiden_->Get())->at(it);
      //if(identk != 10)   continue;
      //std::cout<<identk<<"\n";
      //if(identk != 8)   continue;
      if(identk != 20)   continue;
      double xtk = (xPos_->Get())->at(it);
      double ytk = (yPos_->Get())->at(it);
      tracksPerEvent++;
      //loop over ref hits
      for (unsigned int iref = 0; iref < (sensorId_fp->Get())->size() ; iref++)  {
        int refId = (sensorId_fp->Get())->at(iref);
	//std::cout<<refId<<std::endl;
	//if(refId != 8)    continue;//sensorId 10 is FeI4
        //if(refId != 10)    continue;//sensorId 10 is FeI4
	if(refId != 20)    continue;//sensorId 10 is FeI4
        double xref = (xPos_fp->Get())->at(iref);
        double yref = (yPos_fp->Get())->at(iref);
        //residual fei4
        double xres = xref - xtk;
        double yres = yref - ytk;
        if (std::fabs(xres) < 1. && std::fabs(yres) < 1.) {
          tkhists_->resXfei4->Fill(xres);
          tkhists_->resYfei4->Fill(yres);
          tkhists_->refx_x->Fill(xtk, xres);
          tkhists_->refx_y->Fill(xtk, yres);
          tkhists_->refy_x->Fill(ytk, xres);
          tkhists_->refy_y->Fill(ytk, yres);
        }
        //residual cut at FeI4  
        //if(std::fabs(xres) > 0.25 or std::fabs(yres) > 0.1 )   continue;
	if(std::fabs(xres) > 0.20 or std::fabs(yres) > 0.08 )   continue;
        tkhists_->matchedRefHitMap->Fill(xtk, ytk);
        ref_tracks.push_back((trackNum_->Get())->at(it));
        reftracksPerEvent++;
      }//end loop FeI4
    }//end loop over tracks
    tkhists_->ncleanedTk->Fill(cleanedTk.size());

    //loop over modules
    for(unsigned int imod = 0; imod < nmods_; imod++ ) {
      ModuleHistos& mh = *(mhists_[imod]);
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
      // TOP=32, BOTTOM=33 [Module: 2]
      // row 1 :: FE0
      // row 0 :: FE1
      // for fe0, cbc strip = 1016 - strip
      for (unsigned int nh = 0; nh < (col_->Get())->size() ; nh++)  {
        const int sid =  (iden_->Get())->at(nh);
        const int r = (row_->Get())->at(nh);
	//if (r == 1) continue;
        const int ch = (col_->Get())->at(nh);
        //TOP=30, BOTTOM=31
        //For bottom sensor, mask channel 794 for row==0 [for downstream module only??]
	if (imod == 0) {
	  if( sid == SENSORID::Mod1TOP) {
	    if(r == 1) hitsTopfe0.push_back(ch);
	    else if (r == 0) hitsTopfe1.push_back(ch);
	  }
	  else if( sid == SENSORID::Mod1BOTTOM) { 
	    if(r == 1) hitsBotfe0.push_back(ch);
	    else if (r == 0 and ch != 794) hitsBotfe1.push_back(ch);
	  }
        }
	else if (imod == 1) {
	  if( sid == SENSORID::Mod2TOP) {
	    if(r == 1) hitsTopfe0.push_back(ch);
	    else if (r == 0) hitsTopfe1.push_back(ch);
	  }
	  else if( sid == SENSORID::Mod2BOTTOM) { 
	    if(r == 1) hitsBotfe0.push_back(ch);
	    else if (r == 0) hitsBotfe1.push_back(ch);
	  }
        }
      }//finish loop over hits
      
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


      //test
      for (auto& cls : clsBotfe0) clsBot.push_back(cls);
      for (auto& cls : clsBotfe1) {
	int fstrip = cls.firstStrip() + 1016;
	cluster clust(fstrip, cls.column(), cls.size());
	clsBot.push_back(clust);
      }
      for (auto& cls : clsTopfe0) clsTop.push_back(cls);
      for (auto& cls : clsTopfe1) {
	int fstrip = cls.firstStrip() + 1016;
	cluster clust(fstrip, cls.column(), cls.size());
	clsTop.push_back(clust);
      }
      // correlation
      for(auto& h0: clsBot) {
        for(auto& h1: clsTop)
          mh.clusterPoscorr->Fill(h0.center(), h1.center());
      }

      for (auto& h: clsBot) mh.bottomS.clspos->Fill(h.center());
      for (auto& h: clsTop) mh.topS.clspos->Fill(h.center());
      /////////////////////////////////////////////////			    

      //Plot clus prop
      mh.bottomS.ncls_Fe0->Fill(clsBotfe0.size());
      mh.bottomS.ncls_Fe1->Fill(clsBotfe1.size());
      mh.topS.ncls_Fe0->Fill(clsTopfe0.size());
      mh.topS.ncls_Fe1->Fill(clsTopfe1.size());
      float clsBotfe0_total = 0.0;
      float clsBotfe1_total = 0.0;
      for (auto& h: clsBotfe0) {
	clsBotfe0_total += h.center();
        mh.bottomS.clspos_Fe0->Fill(h.center());
	mh.bottomS.clspos_FeAll->Fill(h.center());
        mh.bottomS.clswidth_Fe0->Fill(h.size());
      }
      clsBotFe0_allevt.push_back(clsBotfe0_total/clsBotfe0.size());
      for (auto& h: clsBotfe1) {
	clsBotfe1_total += h.center();
        mh.bottomS.clspos_Fe1->Fill(h.center());
	mh.bottomS.clspos_FeAll->Fill(h.center());
        mh.bottomS.clswidth_Fe1->Fill(h.size());
      }
      clsBotFe1_allevt.push_back(clsBotfe1_total/clsBotfe1.size());
      for (auto& h: clsTopfe0) {
        mh.topS.clspos_Fe0->Fill(h.center());
	mh.topS.clspos_FeAll->Fill(h.center());
        mh.topS.clswidth_Fe0->Fill(h.size());
      }
      for (auto& h: clsTopfe1) {
        mh.topS.clspos_Fe1->Fill(h.center());
	mh.topS.clspos_FeAll->Fill(h.center());
        mh.topS.clswidth_Fe1->Fill(h.size());
      }
      
      //Correlation histos
      for(auto& h0: clsBotfe0) {
        for(auto& h1: clsTopfe0)
          mh.clusterPoscorrfe0->Fill(h0.center(), h1.center());
      }
      for(auto& h0: clsBotfe1) {
        for(auto& h1: clsTopfe1)
          mh.clusterPoscorrfe1->Fill(h0.center(), h1.center());
      }

      //      std::cout<<"-------Event---------\n";
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
        } else if (fe == 1) cbcstubfe1.emplace_back(stub(strip, fe, bend));
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
        for (auto& h: clsBotfe1) {
          mh.clsvsStubPoscorrfe1->Fill(h.center(), s.center());
        }
      }

      //now do offline stub reco 
      //seeding, matching, stubVec
      stubSimulator(clsBotfe0, clsTopfe0, recostubfe0);//bottom is seeding
      stubSimulator(clsBotfe1, clsTopfe1, recostubfe1);//bottom is seeding

      //plot reco stub prop
      mh.nrecostubsfe0->Fill(recostubfe0.size());
      std::vector<float>recoStubPos;
      for(auto& rs : recostubfe0) {
	float stubPos = (rs.center()-508.5)*0.09;
	recoStubPos.push_back(stubPos);
	mh.recoStubPos -> Fill (stubPos);
        mh.recoStubposfe0->Fill(rs.center());
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
        for (auto& cs: cbcstubfe1) {
          mh.rsvscsPoscorrfe1->Fill( cs.center(), rs.center() );
        }
      }

      //std::cout<<"tdc : "<<tdc<<std::endl;
      //now to Tk to DUT matching
      // If there is ref hit, take the last one
      // if (ref_tracks>-1)
      bool tdcok {false};
      tdcok = true;
      if (imod == 0) tdcok = tdc == 6; // new beamtest
      //if (imod == 0) tdcok = tdc == 2 or tdc == 3; // old beamtest
      if (imod == 1) tdcok = true;
      if(!tdcok)   continue;

      int senID = 0;
      //std::cout<<sensorType_<<"\t"<<senID<<"\n";
      if (sensorType_ == "bottom")   senID = (imod == 0) ? 31 : 33;
      else if (sensorType_ == "top") senID = (imod == 0) ? 30 : 32;
      else {
	std::cout<<"sensortype has not been mentioned correctly! ['bottom' || 'top']\n";
	senID = (imod == 0) ? 31 : 33; //bottom by default
      }
      int nCls = 0;
      bool istop {false};
      bool isbot {false};
      std::vector<double>refTkPos;

      // Loop over reference tracks
      for (unsigned int irt =0; irt < ref_tracks.size(); irt++) {
	int ref_track = ref_tracks[irt]; // get ref track index
	int trackID = -1;
	// Loop over total tracks to check if a track has DUT hit or not
	for (unsigned int it=0; it < (tkiden_->Get())->size(); it++) {
	  // Find the ref track index from the track collection and check if that track 
	  // has a hit at DUT lower or upper sensor
	  if ((trackNum_->Get())->at(it) == ref_track && (tkiden_->Get())->at(it) == senID) {
	    trackID=it;
	    nRefTracks++;
	    if (imod == 0) nRefTracksMod1++;
	    else if (imod == 1) nRefTracksMod2++;
	    break;
	  }
	}
        //Track position on DUT plane
        double txpos = (xPos_->Get())->at(trackID);
        double typos = (yPos_->Get())->at(trackID);

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
        for (unsigned int iref = 0; iref < (sensorId_fp->Get())->size() ; iref++)  {
          int refId = (sensorId_fp->Get())->at(iref);
          if(refId != senID)    continue; //sensorId :: 30,31 || 32,33 is DUT
	  // Get the track fit points on DUT
          double xdut = (xPos_fp->Get())->at(iref);
          double ydut = (yPos_fp->Get())->at(iref);
          unsigned int clsize  = (clustersize_fp->Get())->at(iref);
          //residual DUT
          double xres = xdut - txpos;
          double yres = ydut - typos;
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
	    if (imod == 0) nMatchedTracksMod1++;
	    else if (imod == 1) nMatchedTracksMod2++;
	    mh.nMatchedTk -> Fill(0);
	    nCls++;
            isMatched = true;
	    if (sensorType_ == "top") {
	      mh.topS.clswidth_alltdc->Fill(clsize);
	      if (tdc == 2) mh.topS.clswidth_tdc2->Fill(clsize);
	      else if (tdc == 3) mh.topS.clswidth_tdc3->Fill(clsize);
	      else if (tdc == 4) mh.topS.clswidth_tdc4->Fill(clsize);
	      istop = true;
	    }
	    else if (sensorType_ == "bottom") {
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
	    if (imod == 0) {
	      mh.cbcMatchedTkYPos->Fill(typos);
	      nMatchedcbcStubsMod1++;
	    }
	    else if (imod == 1) nMatchedcbcStubsMod2++;
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
	    if (imod == 0) {
	      nMatchedrecoStubsMod1++;
	      mh.recoMatchedTkYPos->Fill(typos);
	    }
	    else if (imod == 1) nMatchedrecoStubsMod2++;
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
	if (tdc == 2) {
	  mh.bottomS.ncls_tdc2->Fill(nCls);
	}
	else if (tdc == 3) {
	  mh.bottomS.ncls_tdc3->Fill(nCls);
	}
	else if (tdc == 4) {
	  mh.bottomS.ncls_tdc4->Fill(nCls);
	}
      }

      else if (istop) {
	mh.topS.ncls_alltdc->Fill(nCls);//new
	if (tdc == 2) mh.topS.ncls_tdc2->Fill(nCls);//new
	else if (tdc == 3) mh.topS.ncls_tdc3->Fill(nCls);//new
	else if (tdc == 4) mh.topS.ncls_tdc4->Fill(nCls);//new
      }
    }//module loop ends
    //iev++;
  }//event loop while end
 
  raneventloop_ = true;
  std::cout << "Total tracks =" << tracksPerEvent << std::endl; 
  std::cout << "Module [1]\n";
  std::cout << "nRefTracks     = " << nRefTracksMod1 << std::endl;
  std::cout << "nMatchedTracks = " << nMatchedTracksMod1 << std::endl;
  std::cout << "Track matching Efficiency with Clusters [NSA :: 150µ]: " << (float)nMatchedTracksMod1/nRefTracksMod1*100.  << "%" << std::endl;
  std::cout << "nMatchedCBCStubs: "<< nMatchedcbcStubsMod1<<"\t"<< "nMatchedRECOStubs: "<< nMatchedrecoStubsMod1<<"\t"    
	    << "cbcStubEfficiency: "<<(float)nMatchedcbcStubsMod1/nMatchedTracksMod1*100 << "%" << "\t" 
	    << "recoStubEfficiency: "<<(float)nMatchedrecoStubsMod1/nMatchedTracksMod1*100 << "%" <<std::endl;

  std::cout<< "Module [2]\n";
  std::cout << "nRefTracks     = " << nRefTracksMod2 << std::endl;
  std::cout << "nMatchedTracks = " << nMatchedTracksMod2 << std::endl;
  std::cout << "Track matching Efficiency with Clusters [NSA :: 150µ]: " << (float)nMatchedTracksMod2/nRefTracksMod2*100.  << "%" << std::endl;
  std::cout << "nMatchedCBCStubs: "<< nMatchedcbcStubsMod2<<"\t"<< "nMatchedRECOStubs: "<< nMatchedrecoStubsMod2<<"\t"    
	    << "cbcStubEfficiency: "<<(float)nMatchedcbcStubsMod2/nMatchedTracksMod2*100 << "%" << "\t" 
	    << "recoStubEfficiency: "<<(float)nMatchedrecoStubsMod2/nMatchedTracksMod2*100 << "%" <<std::endl;

  /*
  TH2F *hscat = new TH2F("hscat","  ",508,1,1017,508,1,1017);
  for (Int_t i = 0; i < clsBotFe0_allevt.size(); i++) hscat->Fill(clsBotFe0_allevt[i], clsBotFe1_allevt[i]);
  TCanvas *canv = new TCanvas("cst","norm hists",1600,1200);
  canv->cd(1);
  gPad->SetGrid();
  hscat->SetFillColor(kRed);
  hscat->Draw("Box");
  canv->SaveAs("hScat.root");
  */
}//end event loop 

//Function to create clusters from hits
void BTAnalyzer::offlineclusterizer(std::vector<int>& hits, 
				    const unsigned int nCbc, 
				    const unsigned int nStripsPerCBC, 
				    std::vector<cluster>& clusVec ) {
  //std::cout<<"making offline cluster===>>>"<<"\n";
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
  for (unsigned int i = 1; i < hits.size(); i++){
    unsigned int icol    = hits.at(i)/16 < 8 ?  0 : 1;
    //std::cout<<"\t"<<"hit: "<<i<<"\t"<<"strip: "<<hits.at(i)<<"\t"<<"col: "<<icol<<"\n";
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
  //loop over modules
  for(unsigned int i = 0; i < nmods_; i++ ) {  
    ModuleHistos& mh = *(mhists_[i]);
    fout->cd(); 
    mh.writeHistostofile(fout); 
    
  }
  fout->Save();
  fout->Print();
  fout->Close();
  //std::cout<<"bla\n";
}

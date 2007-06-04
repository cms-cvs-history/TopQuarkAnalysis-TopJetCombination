#include <iostream>
#include <cassert>
#include <TROOT.h>
#include <TSystem.h>
#include <Cintex/Cintex.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TH2.h>
#include <TGraph.h>
#include <TF1.h>
#include <TFormula.h>
#include <TStyle.h>
#include <TKey.h>
#include <vector>
#include "FWCore/FWLite/src/AutoLibraryLoader.h"
#include "AnalysisDataFormats/TopObjects/interface/TtSemiEvtSolution.h"
#include "TopQuarkAnalysis/TopTools/interface/LRHelpFunctions.h"

using namespace std;



///////////////////////
// Constants         //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//input files
const  int       nrFiles  	  		= 1;
const  TString   path     	  		= "/beo5/heyninck/CMSSW/src/TopQuarkAnalysis/TopEventProducers/test/TtSemiMuEvents_";

//matching variables
const  bool  	 useSpaceAngle    		= true;
const  double 	 SumAlphaCut  	  		= 0.7;

//loops to be executed
const  bool  	 doJetCombLRObsLoop  		= true;
const  bool  	 doJetCombPurEffLoop  		= true;

//observable histogram variables
const  int      nrJetCombObs  			= 2;
const  int      JetCombObs[nrJetCombObs] 	= {1,3};
const  int   	nrJetCombHistBins    		= 50;
const  double   JetCombObsMin[nrJetCombObs]	= {0,0};
const  double   JetCombObsMax[nrJetCombObs]	= {400,5};

//observable fit functions
const char*     JetCombObsFits[nrJetCombObs] 	= {            //sigmoind
						     "landaun",	//obs0	
						//     "[0]/(1 + 1/exp([1]*([2] - x)))",  //obs1	
						     "[0]/(1 + 1/exp([1]*([2] - x)))"  //obs2	
                                          	  };

//likelihood histogram variables
const  int   	nrJetCombLRtotBins   		= 30;
const  double 	JetCombLRtotMin   		= -10;
const  double 	JetCombLRtotMax      		= -5;
const  char* 	JetCombLRtotFitFunction      	= "[0]/(1 + 1/exp([1]*([2] - x)))";

//output files ps/root
const  TString  JetCombOutfileName   		= "../data/TtSemiLRJetComb.root";
const  TString  JetCombPSfile     		= "../data/TtSemiLRJetComb.ps";

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//
// Global variables
//
LRHelpFunctions *myLRhelper;
void doEventloop(int);
vector<int> obsNrs;
vector<double> obsMin,obsMax;
vector<const char*> obsFits;




//
// Main analysis
//

int main() { 
  gSystem->Load("libFWCoreFWLite");
  AutoLibraryLoader::enable();
  
  
  // define all histograms & fit functions
  //to replace with something more elegant
  for(int j = 0; j < nrJetCombObs; j++){
    obsNrs.push_back(JetCombObs[j]);
    obsMin.push_back(JetCombObsMin[j]);
    obsMax.push_back(JetCombObsMax[j]);
    obsFits.push_back(JetCombObsFits[j]);
  }
  myLRhelper = new LRHelpFunctions(obsNrs, nrJetCombHistBins, obsMin, obsMax, obsFits, 
                                   nrJetCombLRtotBins, JetCombLRtotMin, JetCombLRtotMax, JetCombLRtotFitFunction);
  
  // fill the histograms
  // loop 1: fill signal and background contributions to S and B histograms
  // loop 2: fill calculated LR value for each signal or background contributions
  if(doJetCombLRObsLoop)    { doEventloop(1); myLRhelper -> makeAndFitSoverSplusBHists(); };
  if(! doJetCombLRObsLoop)  myLRhelper -> readObsHistsAndFits(JetCombOutfileName,false);
  if(doJetCombPurEffLoop)   { doEventloop(2); myLRhelper -> makeAndFitPurityHists(); };       
    
  // store histograms and fits in root-file
  myLRhelper -> storeToROOTfile(JetCombOutfileName);
     
  // make some control plots and put them in a .ps file
  myLRhelper -> storeControlPlots(JetCombPSfile);
}





//
// Loop over the events (with the definition of what is considered signal and background)
//

void doEventloop(int loop){ 
  cout<<endl<<endl<<"**** STARTING EVENT LOOP "<<loop<<" ****"<<endl;
  int okEvents = 0, totNrEv = 0;
  for (int nr = 1; nr <= nrFiles; nr++) {
    TString ft = path; ft += nr-1; ft += ".root";
    if (!gSystem->AccessPathName(ft)) {
      TFile *file = TFile::Open(ft);
      TTree * events = dynamic_cast<TTree*>( file->Get( "Events" ) );
      assert( events != 0 );
      TBranch * solsbranch = events->GetBranch( "TtSemiEvtSolutions_solutions__TtEventReco.obj" );
      assert( solsbranch != 0 );
      vector<TtSemiEvtSolution> sols;
      solsbranch->SetAddress( & sols );
    
      //loop over all events in a file 
      for( int ev = 0; ev < events->GetEntries(); ++ ev ) {
        ++totNrEv;
        if((double)((totNrEv*1.)/1000.) == (double) (totNrEv/1000)) cout<< "  Processing event "<< totNrEv<<endl; 
        solsbranch->GetEntry( ev );
        if(sols.size()== 12){
	  double maxLogLRVal = -999.;
	  int    maxLogLRSol = -999;
	  //loop over solutions
	  for(int s=0; s<12; s++){
            // get observable values
	    vector<double> obsVals;
	    for(int j = 0; j < nrJetCombObs; j++){
	      unsigned int o=0;
	      while(sols[s].getLRCorrJetCombVar(o)>0){
	        if(fabs(obsNrs[j]-sols[s].getLRCorrJetCombVar(o))<0.001) obsVals.push_back(sols[s].getLRCorrJetCombVal(o));
	        ++o;
	      }
	    }
	    if(loop==1){
	      // Fill the observables for each jet combination
	      // signal: best MC matching jet combination with a total sumDR of four jets lower than threshold 
	      // background: all other solutions 
	      if(sols[s].getSumDeltaRjp()<SumAlphaCut && sols[s].getMCCorrJetComb()==s) {
	        myLRhelper -> fillToSignalHists(obsVals);
	        ++okEvents;
	      }
	      else
	      {
	        myLRhelper -> fillToBackgroundHists(obsVals);
	      }
            }
	    if(loop==2){
	      double logLR =  myLRhelper -> calcLRval(obsVals);
	      if(logLR>maxLogLRVal) { maxLogLRVal = logLR; maxLogLRSol = s; };
	    }
	  }
	  if(loop==2){
	    if(sols[maxLogLRSol].getSumDeltaRjp()<SumAlphaCut && sols[maxLogLRSol].getMCCorrJetComb()==maxLogLRSol) {
	      myLRhelper -> fillLRSignalHist(maxLogLRVal);
	    }
	    else
	    {
	      myLRhelper -> fillLRBackgroundHist(maxLogLRVal);
	    }
	  }
        }
      }
      file->Close();
    }
    else
    {
      cout<<ft<<" doesn't exist"<<endl;
    }
  }
  if(loop==1){
    cout<<endl<<"***********************  STATISTICS  *************************"<<endl;
    cout<<" Probability that a correct jet combination exists:"<<endl;
    cout<<" (fraction events with ";
    if(useSpaceAngle) cout<<"min SumAngle_jp < ";
    if(!useSpaceAngle) cout<<"min DR_jp < ";
    cout<<SumAlphaCut<<" )"<<endl;
    cout<<endl<<"                 "<<(100.*okEvents)/(1.*totNrEv)<<" %"<<endl;
    cout<<endl<<"******************************************************************"<<endl;
  }
}

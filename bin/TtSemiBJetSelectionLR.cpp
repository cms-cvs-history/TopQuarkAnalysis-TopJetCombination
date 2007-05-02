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
#include "AnalysisDataFormats/TopObjects/interface/TtGenEvent.h"
#include "AnalysisDataFormats/TopObjects/interface/TtSemiEvtSolution.h"
#include "AnalysisDataFormats/TopObjects/interface/BestMatching.h"
#include "TopQuarkAnalysis/TopJetCombination/src/LRHelpFunctions.cc"

using namespace std;
using namespace reco;



//read config file
#include "TopQuarkAnalysis/TopJetCombination/interface/TtJetCombLRsetup.h"



//
// Global variables
//

TH1F     * hBSelObs[nrBSelObs][3], * hBSelLRtotS, * hBSelLRtotB, * hBSelPurity, * hAngleDiff; 
TGraph   * hBSelEffvsPur;
TH2F     * corBSelObs[nrBSelObs][nrBSelObs];
TF1      * bSelFits[nrBSelObs], * bSelFits2[nrBSelObs], * fBSelPurity, combBTagFit; 
TCanvas  * myBSelC[nrBSelObs];
TFile    * bSelOutfile;

void fillObservableValues(TtGenEvent,vector<TtSemiEvtSolution>,vector<double>);
void fillLRtoHistos(TtGenEvent,vector<TtSemiEvtSolution>,vector<double>);
double LRtotFMin;
double LRtotFMax;







//
// Main analysis
//

int main() {
  gSystem->Load( "libFWCoreFWLite" );
  AutoLibraryLoader::enable();
  
  if(doBSelLRObsLoop){
    combBTagFit = readBTagLR("../data/TtSemiBProbLR.root");
  
    //define histograms
    hAngleDiff = new TH1F("hAngleDiff","hAngleDiff",60,-1,4);
    if(useSpaceAngle)  hAngleDiff->GetXaxis()->SetTitle("min. sum Alpha_jp");
    if(!useSpaceAngle) hAngleDiff->GetXaxis()->SetTitle("min. sum DR_jp");
    for (int j = 0; j < nrBSelObs; j++) {
      for (int i = 0; i < 2; i++) {
        hBSelObs[j][i] = new TH1F("dummy", "", nrBSelHistBins, (double) bSelObsMin[j], (double) bSelObsMax[j]);
        TString ht = "hBSelObs"; ht += bSelObs[j]; ht += "_"; 
	if(i==0) ht += "S"; if(i==1) ht += "B"; 
        hBSelObs[j][i] -> SetName(ht);
        hBSelObs[j][i] -> GetXaxis() -> SetTitle(ht);
      }
      for (int i = j; i < nrBSelObs; i++) {
        corBSelObs[j][i] = new TH2F("dummy", "Correlation histogram", nrBSelHistBins, (double) bSelObsMin[j], (double) bSelObsMax[j], nrBSelHistBins, (double) bSelObsMin[i],(double)  bSelObsMax[i]);
        TString ht = "corBSelObs"; ht += bSelObs[j]; ht += "_Obs"; ht += bSelObs[i]; 
        corBSelObs[j][i] -> SetName(ht);
        TString t1 = "Obs"; t1 += bSelObs[j]; corBSelObs[j][i] -> GetXaxis() -> SetTitle(t1);
        TString t2 = "Obs"; t2 += bSelObs[i]; corBSelObs[j][i] -> GetYaxis() -> SetTitle(t2);
      }
    }
  
    //loop over events to fill observables
    int okEvents = 0;
    for (int nr = 1; nr <= nrFiles; nr++) {
      TString ft = path; 
      ft += nr; ft += ".root";
      if (!gSystem->AccessPathName(ft)) {
        TFile *file = TFile::Open(ft);
        TTree * events = dynamic_cast<TTree*>( file->Get( "Events" ) );
        assert( events != 0 );
        TBranch * solsbranch = events->GetBranch( "TtSemiEvtSolutions_fitsolutions_iterativeCone5CaloJets_TTEVENT.obj" );
        assert( solsbranch != 0 );
        TBranch * partonbranch = events->GetBranch( "TtGenEvent_genEvt__TTEVENT.obj" );
        assert( partonbranch != 0 );
        vector<TtSemiEvtSolution> sols;
        solsbranch->SetAddress( & sols );
        TtGenEvent genEvt;
        partonbranch->SetAddress( & genEvt );

        int nev = events->GetEntries();
    
        for( int ev = 0; ev < nev; ++ ev ) {
          solsbranch->GetEntry( ev );
          partonbranch->GetEntry( ev );
          if(genEvt.decay()==1){
	    if(sols.size()== 12)
	    {
	      vector<double> bestMatch = BestMatch(genEvt, sols, useSpaceAngle);
	      hAngleDiff->Fill(bestMatch[2]);
	      if(bestMatch[2]<SumAlphaCut) {
	        fillObservableValues(genEvt,sols,bestMatch); //SemiLep decay
		++okEvents;
              }
	    }
	    else
	    {
	      hAngleDiff->Fill(-1.);
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
    
    //get observable fitting functions
    vector<TF1> fits = getBSelObsFitFunctions();
    for (int i = 0; i < nrBSelObs; i++) {
      bSelFits[i] = new TF1(fits[bSelObs[i]]);
      bSelFits[i] -> SetRange(bSelObsMin[i],bSelObsMax[i]);
    }
    
    //make Signal over Signal + Background Plots
    for (int i = 0; i < nrBSelObs; i++) {
      hBSelObs[i][2] = new TH1F(SoverB(hBSelObs[i][0],hBSelObs[i][1],bSelObs[i]));
      hBSelObs[i][2] -> Fit(bSelFits[i],"RQ");
    }     
    
    
    //output file
    bSelOutfile = new TFile(BSeloutfileName, "RECREATE");
    bSelOutfile->cd();
    hAngleDiff -> Write();
    for (int j = 0; j < nrBSelObs; j++) {
      for (int i = 0; i < 3; i++) {
        hBSelObs[j][i] -> Write();
      }
      for (int i = j; i < nrBSelObs; i++) {
        corBSelObs[j][i] -> Write();
      }
    }
    bSelOutfile->Close();
    delete bSelOutfile;
    
    cout<<endl<<"******************************************************************"<<endl;
    cout<<" Probability that a correct jet combination exists:"<<endl;
    cout<<" (fraction events with ";
    if(useSpaceAngle) cout<<"min SumAngle_jp < ";
    if(!useSpaceAngle) cout<<"min DR_jp < ";
    cout<<SumAlphaCut<<" )"<<endl;
    cout<<endl<<"                 "<<(100.*okEvents)/(1.*hAngleDiff->GetEntries())<<" %"<<endl;
    cout<<endl<<"******************************************************************"<<endl;
    
    cout<<"Frac given by plot method = "<<getPExistingTrueComb(BSeloutfileName)<<endl;
  }
  
  
  
  
  if(doBSelPurEffLoop){  
  
    vector<TF1> fits = getBSelObsFits(BSeloutfileName);
    for (int i = 0; i < nrBSelObs; i++) {
      bSelFits2[i] = new TF1(fits[i]);
    }

    //define some histograms
    hBSelLRtotS = new TH1F("hBSelLRtotS","",nrBSelLRtotBins,BSelLRtotMin,BSelLRtotMax);
    hBSelLRtotS -> GetXaxis() -> SetTitle("log combined LR for signal");
    hBSelLRtotB = new TH1F("hBSelLRtotB","",nrBSelLRtotBins,BSelLRtotMin,BSelLRtotMax);
    hBSelLRtotB -> GetXaxis() -> SetTitle("log combined LR for background");
    fBSelPurity = new TF1("fBSelPurity","[0]-[1]/(1 + 1/exp([2]*(x-[3])))");
    fBSelPurity->SetParameters(0.4,-0.6,2,-1);
    
    //make Purity Plot
    LRtotFMin = 1000.;
    LRtotFMax = -1000.;
    //loop over events to fill observables
    for (int nr = 1; nr <= nrFiles; nr++) {
      TString ft = path; 
      ft += nr; ft += ".root";
      if (!gSystem->AccessPathName(ft)) {
        TFile *file = TFile::Open(ft);
        TTree * events = dynamic_cast<TTree*>( file->Get( "Events" ) );
        assert( events != 0 );
        TBranch * solsbranch = events->GetBranch( "TtSemiEvtSolutions_fitsolutions_iterativeCone5CaloJets_TTEVENT.obj" );
        assert( solsbranch != 0 );
        TBranch * partonbranch = events->GetBranch( "TtGenEvent_genEvt__TTEVENT.obj" );
        assert( partonbranch != 0 );
        vector<TtSemiEvtSolution> sols;
        solsbranch->SetAddress( & sols );
        TtGenEvent genEvt;
        partonbranch->SetAddress( & genEvt );

        int nev = events->GetEntries();
    
        for( int ev = 0; ev < nev; ++ ev ) {
          solsbranch->GetEntry( ev );
          partonbranch->GetEntry( ev );
          if(genEvt.decay()==1 && sols.size()== 12)
	  {
	      vector<double> bestMatch = BestMatch(genEvt, sols, useSpaceAngle);
	      if(bestMatch[2]<SumAlphaCut) fillLRtoHistos(genEvt,sols,bestMatch); //SemiLep decay
          }
        }
        file->Close();
      }
      else
      {
        cout<<ft<<" doesn't exist"<<endl;
      }
    } 
      
      
    cout<<"Range combined LR = " <<LRtotFMin<<" -> "<<LRtotFMax<<endl;
    
    //Make Purity plots
    hBSelPurity   = new TH1F(makePurityPlot(hBSelLRtotS, hBSelLRtotB)); 
    hBSelPurity   -> Fit(fBSelPurity,"Q");  
    hBSelEffvsPur = new TGraph(makeEffVsPurGraph(hBSelLRtotS,fBSelPurity));
    
    //re-open output file & save extra plots
    bSelOutfile = new TFile(BSeloutfileName, "UPDATE");
    bSelOutfile->cd();
    hBSelLRtotS -> Write();
    hBSelLRtotB -> Write();
    hBSelPurity -> Write();
    hBSelEffvsPur -> Write();
    bSelOutfile->Close();

  }
}












//////////////////////////////////////////////////////////////////////
//
// Fill LR Observable values for 1 event
//////////////////////////////////////////////////////////////////////

void fillObservableValues(TtGenEvent genEvt,vector<TtSemiEvtSolution> sols, vector<double> match) {

  double sumCorrectBjetsEnergy = sols[(int)match[1]].getHadb().energy()+sols[(int)match[1]].getLepb().energy();
  for(unsigned int s=0; s<sols.size(); s++){
    vector<double> bSelObsVal = getBSelObsValues(&sols[s],&combBTagFit);
    double sumBjetsEnergy = sols[s].getHadb().energy()+sols[s].getLepb().energy();
    if(fabs(sumCorrectBjetsEnergy-sumBjetsEnergy)<0.001)
    {
      for (int j = 0; j < nrBSelObs; j++) {
        hBSelObs[j][0]->Fill(bSelObsVal[bSelObs[j]]);
        for (int k = j; k < nrBSelObs; k++) {
	  corBSelObs[j][k] -> Fill(bSelObsVal[bSelObs[j]],bSelObsVal[bSelObs[k]]);
	}
      }
    }
    else
    {
      for (int j = 0; j < nrBSelObs; j++) {
        hBSelObs[j][1]->Fill(bSelObsVal[bSelObs[j]]);
      }
    }
  }   
}













//////////////////////////////////////////////////////////////////////////////////////////////
//
// Fill combined LR value for 1 event
//////////////////////////////////////////////////////////////////////////////////////////////

void fillLRtoHistos(TtGenEvent genEvt,vector<TtSemiEvtSolution> sols, vector<double> match){

  double sumCorrectBjetsEnergy = sols[(int)match[1]].getHadb().energy()+sols[(int)match[1]].getLepb().energy();
  
  for(unsigned int s=0; s<sols.size(); s++){
    vector<double> bSelObsVal = getBSelObsValues(&sols[s],&combBTagFit);
    
    //calculate logLR value for this jet combination solution
    double logLR = 0;
    for (int j = 0; j < nrBSelObs; j++) {
      logLR += log(bSelFits2[j]->Eval(bSelObsVal[bSelObs[j]]));
    }
    
    //Fill logLR value
    double sumBjetsEnergy = sols[s].getHadb().energy()+sols[s].getLepb().energy();
    if(fabs(sumCorrectBjetsEnergy-sumBjetsEnergy)<0.001)
    {
      hBSelLRtotS -> Fill(logLR);
      LRtotFMin = fmin(LRtotFMin,logLR);
      LRtotFMax = fmax(LRtotFMax,logLR);
    }
    else
    {
      hBSelLRtotB -> Fill(logLR);
      LRtotFMin = fmin(LRtotFMin,logLR);
      LRtotFMax = fmax(LRtotFMax,logLR);
    }
  }
  
}




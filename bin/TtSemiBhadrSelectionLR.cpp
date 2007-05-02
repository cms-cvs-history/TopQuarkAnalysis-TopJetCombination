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
#include <TKey.h>
#include <TFormula.h>
#include <TStyle.h>
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

TH1F     * hBhadrObs[nrBhadrObs][3], * hBhadrLRtotS, * hBhadrLRtotB, * hBhadrPurity; 
TGraph   * hBhadrEffvsPur;
TH2F     * corBhadrObs[nrBhadrObs][nrBhadrObs];
TF1      * bHadrFits[nrBhadrObs], * bHadrFits2[nrBhadrObs], * fBhadrPurity, combBTagFit; 
TCanvas  * myBhadrC[nrBhadrObs];
TFile    * bHadrOutfile;

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
  
  if(doBhadrLRObsLoop){
    combBTagFit = readBTagLR("../data/TtSemiBProbLR.root");
  
    //define histograms
    for (int j = 0; j < nrBhadrObs; j++) {
      for (int i = 0; i < 2; i++) {
        hBhadrObs[j][i] = new TH1F("dummy", "", nrBhadrHistBins, (double) bHadrObsMin[j], (double) bHadrObsMax[j]);
        TString ht = "hBhadrObs"; ht += bHadrObs[j]; ht += "_"; 
	if(i==0) ht += "S"; if(i==1) ht += "B"; 
        hBhadrObs[j][i] -> SetName(ht);
        hBhadrObs[j][i] -> GetXaxis() -> SetTitle(ht);
      }
      for (int i = j; i < nrBhadrObs; i++) {
        corBhadrObs[j][i] = new TH2F("dummy", "Correlation histogram", nrBhadrHistBins, (double) bHadrObsMin[j], (double) bHadrObsMax[j], nrBhadrHistBins, (double) bHadrObsMin[i],(double)  bHadrObsMax[i]);
        TString ht = "corBhadrObs"; ht += bHadrObs[j]; ht += "_Obs"; ht += bHadrObs[i]; 
        corBhadrObs[j][i] -> SetName(ht);
        TString t1 = "Obs"; t1 += bHadrObs[j]; corBhadrObs[j][i] -> GetXaxis() -> SetTitle(t1);
        TString t2 = "Obs"; t2 += bHadrObs[i]; corBhadrObs[j][i] -> GetYaxis() -> SetTitle(t2);
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
	      if(bestMatch[2]<SumAlphaCut) {
	        fillObservableValues(genEvt,sols,bestMatch); //SemiLep decay
		++okEvents;
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
    
    //get observable fitting functions
    vector<TF1> fits = getBhadrObsFitFunctions();
    for (int i = 0; i < nrBhadrObs; i++) {
      bHadrFits[i] = new TF1(fits[bHadrObs[i]]);
      bHadrFits[i] -> SetRange(bHadrObsMin[i],bHadrObsMax[i]);
    }
    
    //make Signal over Signal + Background Plots
    for (int i = 0; i < nrBhadrObs; i++) {
      hBhadrObs[i][2] = new TH1F(SoverB(hBhadrObs[i][0],hBhadrObs[i][1],bHadrObs[i]));
      hBhadrObs[i][2] -> Fit(bHadrFits[i],"RQ");
    }     
    
    
    //output file
    bHadrOutfile = new TFile(BhadroutfileName, "RECREATE");
    bHadrOutfile->cd();
    for (int j = 0; j < nrBhadrObs; j++) {
      for (int i = 0; i < 3; i++) {
        hBhadrObs[j][i] -> Write();
      }
      for (int i = j; i < nrBhadrObs; i++) {
        corBhadrObs[j][i] -> Write();
      }
    }
    bHadrOutfile->Close();
    delete bHadrOutfile;
    
  }
  
  
  
  
  if(doBhadrPurEffLoop){  
  
    vector<TF1> fits = getBhadrObsFits(BhadroutfileName);
    for (int i = 0; i < nrBhadrObs; i++) {
      bHadrFits2[i] = new TF1(fits[i]);
    }

    //define some histograms
    hBhadrLRtotS = new TH1F("hBhadrLRtotS","",nrBhadrLRtotBins,BhadrLRtotMin,BhadrLRtotMax);
    hBhadrLRtotS -> GetXaxis() -> SetTitle("log combined LR for signal");
    hBhadrLRtotB = new TH1F("hBhadrLRtotB","",nrBhadrLRtotBins,BhadrLRtotMin,BhadrLRtotMax);
    hBhadrLRtotB -> GetXaxis() -> SetTitle("log combined LR for background");
    fBhadrPurity = new TF1("fBhadrPurity","[0]-[1]/(1 + 1/exp([2]*(x-[3])))");
    fBhadrPurity->SetParameters(0.4,-0.6,2,-1);
    
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
    hBhadrPurity   = new TH1F(makePurityPlot(hBhadrLRtotS, hBhadrLRtotB)); 
    hBhadrPurity   -> Fit(fBhadrPurity,"Q");  
    hBhadrEffvsPur = new TGraph(makeEffVsPurGraph(hBhadrLRtotS,fBhadrPurity));
    
    //re-open output file & save extra plots
    bHadrOutfile = new TFile(BhadroutfileName, "UPDATE");
    bHadrOutfile->cd();
    hBhadrLRtotS -> Write();
    hBhadrLRtotB -> Write();
    hBhadrPurity -> Write();
    hBhadrEffvsPur -> Write();
    bHadrOutfile->Close();

  }
}












//////////////////////////////////////////////////////////////////////
//
// Fill LR Observable values for 1 event
//////////////////////////////////////////////////////////////////////

void fillObservableValues(TtGenEvent genEvt,vector<TtSemiEvtSolution> sols, vector<double> match) {

  double sumCorrectBjetsEnergy = sols[(int)match[1]].getHadb().energy()+sols[(int)match[1]].getLepb().energy();
  for(unsigned int s=0; s<sols.size(); s++){
    vector<double> bHadrObsVal = getBhadrObsValues(&sols[s],&combBTagFit);
    double sumBjetsEnergy = sols[s].getHadb().energy()+sols[s].getLepb().energy();
    if(fabs(sumCorrectBjetsEnergy-sumBjetsEnergy)<0.001)
    {
      if((unsigned int)match[1] == s){
        for (int j = 0; j < nrBhadrObs; j++) {
          hBhadrObs[j][0]->Fill(bHadrObsVal[bHadrObs[j]]);
          for (int k = j; k < nrBhadrObs; k++) {
	    corBhadrObs[j][k] -> Fill(bHadrObsVal[bHadrObs[j]],bHadrObsVal[bHadrObs[k]]);
	  }
        }
      }
      else
      {
        for (int j = 0; j < nrBhadrObs; j++) {
          hBhadrObs[j][1]->Fill(bHadrObsVal[bHadrObs[j]]);
          for (int k = j; k < nrBhadrObs; k++) {
	    corBhadrObs[j][k] -> Fill(bHadrObsVal[bHadrObs[j]],bHadrObsVal[bHadrObs[k]]);
	  }
        }
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
    vector<double> bHadrObsVal = getBhadrObsValues(&sols[s],&combBTagFit);
    
    //calculate logLR value for this jet combination solution
    double logLR = 0;
    for (int j = 0; j < nrBhadrObs; j++) {
      logLR += log(bHadrFits2[j]->Eval(bHadrObsVal[bHadrObs[j]]));
    }
    
    //Fill logLR value
    double sumBjetsEnergy = sols[s].getHadb().energy()+sols[s].getLepb().energy();
    if(fabs(sumCorrectBjetsEnergy-sumBjetsEnergy)<0.001)
    {
      if((unsigned int)match[1] == s){
        hBhadrLRtotS -> Fill(logLR);
        LRtotFMin = fmin(LRtotFMin,logLR);
        LRtotFMax = fmax(LRtotFMax,logLR);
      }
      else
      {
        hBhadrLRtotB -> Fill(logLR);
        LRtotFMin = fmin(LRtotFMin,logLR);
        LRtotFMax = fmax(LRtotFMax,logLR);
      }
    }
  }
  
}




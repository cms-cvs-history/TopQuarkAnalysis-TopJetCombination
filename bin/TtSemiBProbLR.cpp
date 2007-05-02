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
#include <TKey.h>
#include <TF1.h>
#include <TFormula.h>
#include <TStyle.h>
#include <vector>
#include "FWCore/FWLite/src/AutoLibraryLoader.h"
#include "AnalysisDataFormats/TopObjects/interface/TtGenEvent.h"
#include "AnalysisDataFormats/TopObjects/interface/TtSemiEvtSolution.h"

using namespace std;
using namespace reco;

//config file
#include "TopQuarkAnalysis/TopJetCombination/interface/TtJetCombLRsetup.h"







TH1F * obsHist[4];
TFile * outfile;
void defineParticleHistos();
void drawPartHistos(TCanvas * myCanvas, TH1F * myHists[4]);
void fillSoverSB(TH1F * myHists[4]);
void doParticleOutput();
void AnalyseParticle(TtGenEvent genEvt, TtSemiEvtSolution sol);

//
// Main code
//
int main() {
  gSystem->Load( "libFWCoreFWLite" );
  AutoLibraryLoader::enable();

  outfile = new TFile("../data/TtSemiBProbLR.root", "RECREATE");;
  defineParticleHistos();

  for (int nr = 1; nr <= nrFiles; nr++) {
    TString ft = path; 
    ft += nr; ft += ".root";
    if (!gSystem->AccessPathName(ft)) {
      TFile *file = TFile::Open(ft);
      TTree * events = dynamic_cast<TTree*>( file->Get( "Events" ) );
      assert( events != 0 );
      TBranch * solsbranch = events->GetBranch( "TtSemiEvtSolutions_calsolutions_iterativeCone5CaloJets_TTEVENT.obj" );
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
        if(genEvt.decay()==1 && sols.size()== 12) AnalyseParticle(genEvt,sols[0]); //SemiLep decay
      }
      file->Close();
    }
    else
    {
      cout<<ft<<" doesn't exist"<<endl;
    }
  }
  
  gROOT->SetBatch();
  gROOT->SetStyle("Plain");  
  fillSoverSB(obsHist);
  doParticleOutput();

  outfile->cd();
  outfile->Write();
  outfile->Close();

  return 0;
}


//
// create and define the histograms
//
void defineParticleHistos() {
  for (int i = 0; i < 4; i++) {
      obsHist[i] = new TH1F("dummy", "Observable histogram", 250, -3, 10);
      TString title = "obsHist"; title += i;
      obsHist[i]->SetName(title);
  }
}

	
//
// Analysis for 1 event
//
void AnalyseParticle(TtGenEvent genEvt, TtSemiEvtSolution sol) {
 
  vector<TtJet> jets;
  jets.push_back(sol.getHadp());
  jets.push_back(sol.getHadq());
  jets.push_back(sol.getHadb());
  jets.push_back(sol.getLepb());

  for(unsigned int j = 0; j < jets.size(); j++) {
    bool isol = true;
    for (unsigned int k = 0; k < jets.size(); k++) {
        double phidiff = fabs(jets[j].phi()-jets[k].phi());
        if(phidiff>3.14159) phidiff = phidiff -3.14159;
        double deltaR = sqrt(pow(jets[j].eta()-jets[k].eta(),2)+pow(phidiff,2));
        if (j != k && deltaR < .3) isol = false;
    }
    if (isol) {
      bool fromGenb = false;
      double phidiff = fabs(jets[j].phi()-genEvt.particles()[2].momentum().phi());
      if(phidiff>3.14159) phidiff = phidiff -3.14159;
      double deltaR = sqrt(pow(jets[j].eta()-genEvt.particles()[2].momentum().eta(),2)+pow(phidiff,2));
      if(deltaR <0.3) fromGenb = true;
      phidiff = fabs(jets[j].phi()-genEvt.particles()[3].momentum().phi());
      if(phidiff>3.14159) phidiff = phidiff -3.14159;
      deltaR = sqrt(pow(jets[j].eta()-genEvt.particles()[3].momentum().eta(),2)+pow(phidiff,2));
      if(deltaR <0.3) fromGenb = true;
      if (fromGenb) {
        obsHist[0]->Fill(jets[j].getBdiscriminant());
      } else {
        obsHist[2]->Fill(jets[j].getBdiscriminant());
      }
      //cout << "isolated jet" << endl;
    } else {
      //cout << "non isolated jet" << endl;
      obsHist[1]->Fill(-10);
    }
  }
}



//
// Function to draw the particle based histograms
//
void drawPartHistos(TCanvas * myCanvas, TH1F * myHists[4]) {
  myCanvas->Divide(2,2);
  myCanvas->cd(1); myHists[0]->Draw();
  myCanvas->cd(3); myHists[2]->Draw();
  myCanvas->cd(4); myHists[3]->Draw();
}



//
// Function to draw the Signal over Signal+Background plot
//
void fillSoverSB(TH1F * myHists[4]) {
  for (int k = 1; k <= 250; k++) {
    if (myHists[0]->GetBinContent(k) == 0 && myHists[2]->GetBinContent(k) == 0) {
      myHists[3]->SetBinContent(k, 0);
      myHists[3]->SetBinError(k, 0);
    } else {
      myHists[3]->SetBinContent(k, myHists[0]->GetBinContent(k) / (myHists[2]->GetBinContent(k) + myHists[0]->GetBinContent(k)));
      myHists[3]->SetBinError(k, sqrt(myHists[0]->GetBinContent(k) * myHists[2]->GetBinContent(k) / pow(myHists[0]->GetBinContent(k) + myHists[2]->GetBinContent(k), 3)));
    }
  }
}



//
// Particle Output
//
void doParticleOutput() {

  // define some nice fit functions
  TFormula gauss("gauss", "gaus");
  TFormula symgauss("symgauss", "[0]*(exp(-0.5*(x/[1])**2))");
  TFormula dblgauss("dblgauss", "[0]*(exp(-0.5*((x-[1])/[2])**2)+exp(-0.5*((x+[3])/[4])**2))");
  TFormula symdblgauss("symdblgauss", "[0]*(exp(-0.5*((x-[1])/[2])**2)+exp(-0.5*((x+[1])/[2])**2))");
  TFormula sigm("sigm", "[0]/(1 + 1/exp([1]*([2] - x)))");
  TFormula sigmc("sigmc", "[0]/(1 + 1/exp([1]*([2] - x)))+[3]");
  TFormula dblsigm("dblsigm", "[0]/(1 + 1/exp([1]**2*([2] - x)))/(1 + 1/exp([3]**2*(x - [4])))");
  TFormula symdblsigm("symdblsigm", "[0]/(1 + 1/exp([1]**2*([2] - x)))/(1 + 1/exp([1]**2*([2] + x)))");

  TCanvas *myC = new TCanvas("myC", "Observable canvas");
  drawPartHistos(myC, obsHist); 
  TF1 * obsFit = new TF1("obsFit", "pol2+sigm");
  obsFit->SetParameters(.2, .1, -.007, .5, -3, -0.1);
  obsHist[3]->Fit(obsFit);

}

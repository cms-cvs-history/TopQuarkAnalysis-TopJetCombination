//
// Author:  Jan Heyninck
// Created: Tue Apr  3 17:33:23 PDT 2007
//
// $Id: TtSemiLRJetCombObservables.cc,v 1.6.8.3 2008/04/14 22:35:42 rwolf Exp $
//
#include "TopQuarkAnalysis/TopJetCombination/interface/TtSemiLRJetCombObservables.h"
#include "PhysicsTools/Utilities/interface/DeltaR.h"
// #include "PhysicsTools/Utilities/interface/DeltaPhi.h"
// #include "PhysicsTools/Utilities/interface/DeltaTheta.h"
#include "AnalysisDataFormats/TopObjects/interface/TopJet.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Common/interface/Handle.h"
#include "AnalysisDataFormats/TopObjects/interface/TtGenEvent.h"

using namespace reco;
using namespace std;
using namespace math;
using namespace edm;

// constructor with path; default should not be used
TtSemiLRJetCombObservables::TtSemiLRJetCombObservables() {}


// destructor
TtSemiLRJetCombObservables::~TtSemiLRJetCombObservables() {}

std::vector< TtSemiLRJetCombObservables::IntBoolPair >
TtSemiLRJetCombObservables::operator() (TtSemiEvtSolution &solution, const edm::Event & iEvent, bool matchOnly)
{ 
  bool debug=false;
  
  evtselectVarVal.clear();
  evtselectVarMatch.clear();
 
  // Check whether the objects are matched:
  bool matchHadt = false;
  bool matchLept = false;
  bool matchHadW = false;
  bool matchLepW = false;
  bool matchHadb = false;
  bool matchLepb = false;
  bool matchHadp = false;
  bool matchHadq = false;  
  bool matchHadpq = false;
  bool matchHadqp = false;
  bool matchLepl = false;
  bool matchLepn = false;
  
  if(debug) cout << "== start matching the objects " << endl;
  
  try {
    
    if(debug) cout << "== start trying " << endl;
    double drLepl=0, drLepn=0, drHadb=0, drLepb=0, drHadp=0, drHadq=0, drHadpq=0, drHadqp=0, drHadt=0, drLept=0, drLepW=0, drHadW=0;
    
    Handle<TtGenEvent> genEvent;
    iEvent.getByLabel ("genEvt",genEvent);
    if(debug) cout << "== found genEvent " << endl;
    
    if (genEvent->isSemiLeptonic() && genEvent->numberOfBQuarks() == 2) {

      //if(debug) cout << "== genEvent->quarkFromAntiTop() " << genEvent->quarkFromAntiTop()->pt() << endl;
      if(debug) cout << "== genEvent->isSemiLeptonic() && genEvent->numberOfBQuarks() == 2 " << endl;
      if(debug) cout << "== solution.getDecay() == " <<solution.getDecay()<< endl;     
      if(debug) cout << "== solution.getRecLepm().pt() = " <<solution.getRecLepm().pt()  << endl;
      //if(debug) if(solution.getGenLepl() == 0) cout << "solution.getGenLepl() == NULL" << endl;
      if(debug) cout << "== *(solution.getGenLept())" << solution.getGenLept()->pt() << endl;
      if(debug) cout << "== *(solution.getGenLepl())" << solution.getGenLepl()->pt() << endl;
      // cout << "Semilepton:\n";
      // Match the lepton by deltaR
      if (solution.getDecay() == "muon")     drLepl = DeltaR<reco::Particle>()(solution.getRecLepm(), *(solution.getGenLepl()));
      if(debug) cout << "== matching lepton " << endl;
      if (solution.getDecay() == "electron") drLepl = DeltaR<reco::Particle>()(solution.getRecLepe(), *(solution.getGenLepl()));
      matchLepl = (drLepl < 0.3);
      
      if(debug) cout << "== lepton is matched " << endl;
      // Match the neutrino by deltaR
      drLepn = DeltaR<reco::Particle>()(solution.getRecLepn(), *(solution.getGenLepn()));
      matchLepn = (drLepn < 0.3);

      // Match the hadronic b by deltaR
      drHadb = DeltaR<reco::Particle>()(solution.getRecHadb(), *(solution.getGenHadb()));
      matchHadb = (drHadb < 0.3);

      // Match the hadronicleptonic b by deltaR
      drLepb = DeltaR<reco::Particle>()(solution.getRecLepb(), *(solution.getGenLepb()));
      matchLepb = (drLepb < 0.3);

      // Match the hadronic p by deltaR
      drHadp = DeltaR<reco::Particle>()(solution.getRecHadp(), *(solution.getGenHadp()));
      matchHadp = (drHadp < 0.3);
      
      // Match the hadronic pq by deltaR
      drHadpq = DeltaR<reco::Particle>()(solution.getRecHadp(), *(solution.getGenHadq()));
      matchHadpq = (drHadpq < 0.3);
    
      // Match the hadronic q by deltaR
      drHadq = DeltaR<reco::Particle>()(solution.getRecHadq(), *(solution.getGenHadq()));
      matchHadq = (drHadq < 0.3);      

      // Match the hadronic qp by deltaR
      drHadqp = DeltaR<reco::Particle>()(solution.getRecHadq(), *(solution.getGenHadp()));
      matchHadqp = (drHadqp < 0.3);  

      // Match the hadronic W by deltaR
      drHadW = DeltaR<reco::Particle>()(solution.getRecHadW(), *(solution.getGenHadW()));
      matchHadW = (drHadW < 0.3);    

      // Match the leptonic W by deltaR
      drLepW = DeltaR<reco::Particle>()(solution.getRecLepW(), *(solution.getGenLepW()));
      matchLepW = (drLepW < 0.3);  
     
      // Match the hadronic t by deltaR
      drHadt = DeltaR<reco::Particle>()(solution.getRecHadt(), *(solution.getGenHadt()));
      matchHadt = (drHadt < 0.3);    

      // Match the leptonic t by deltaR
      drLept = DeltaR<reco::Particle>()(solution.getRecLept(), *(solution.getGenLept()));
      matchLept = (drLept < 0.3);   
    }
  } catch (...){cout << "Exception\n";}

  if(debug) cout << "== objects matched" <<endl;

  Handle<vector<TopJet> > jets;
  iEvent.getByLabel(jetSource_, jets);

  if(debug) cout << "== start calculating observables" << endl;

  //obs1 : pt(had top) 
  double AverageTop =((solution.getHadb().p4()+solution.getHadq().p4()+solution.getHadp().p4()).pt()+(solution.getLepb().p4()+solution.getHadq().p4()+solution.getHadp().p4()).pt()+(solution.getHadb().p4()+solution.getLepb().p4()+solution.getHadp().p4()).pt()+(solution.getHadb().p4()+solution.getHadq().p4()+solution.getLepb().p4()).pt())/4.;
  double Obs1 = ((solution.getHadb().p4()+solution.getHadq().p4()+solution.getHadp().p4()).pt())/AverageTop;
  evtselectVarVal.push_back(IntDblPair(1,Obs1));
  evtselectVarMatch.push_back(IntBoolPair(1, ((matchHadq&&matchHadp)||(matchHadpq&&matchHadqp))&&matchHadb)); 
 
  //obs2 : (pt_b1 + pt_b2)/(sum jetpt)
  double Obs2 = (solution.getHadb().pt()+solution.getLepb().pt())/(solution.getHadp().pt()+solution.getHadq().pt());
  evtselectVarVal.push_back(IntDblPair(2,Obs2));
  evtselectVarMatch.push_back(IntBoolPair(2,((matchHadp&&matchHadq)||(matchHadpq&&matchHadqp))&&matchHadb&&matchLepb)); 

  //obs3: delta R between lep b and lepton 
  double Obs3 = -10;
  if (solution.getDecay() == "muon")     Obs3 = ROOT::Math::VectorUtil::DeltaR( solution.getLepb().p4(),solution.getRecLepm().p4() );
  if (solution.getDecay() == "electron") Obs3 = ROOT::Math::VectorUtil::DeltaR( solution.getLepb().p4(),solution.getRecLepe().p4() );
  evtselectVarVal.push_back(IntDblPair(3,Obs3));
  evtselectVarMatch.push_back(IntBoolPair(3,matchLepb&&matchLepl)); 
  
   //obs4 : del R ( had b, had W)
  double Obs4 = ROOT::Math::VectorUtil::DeltaR( solution.getHadb().p4(), solution.getHadq().p4()+solution.getHadp().p4() );
  evtselectVarVal.push_back(IntDblPair(4,Obs4));  
  evtselectVarMatch.push_back(IntBoolPair(4,matchHadb&&((matchHadp&&matchHadp)||(matchHadpq&&matchHadqp)))); 
 
  //obs5 : del R between light quarkssolution.getHadp().phi(
  double Obs5 = ROOT::Math::VectorUtil::DeltaR( solution.getHadq().p4(),solution.getHadp().p4() );
  evtselectVarVal.push_back(IntDblPair(5,Obs5)); 
  evtselectVarMatch.push_back(IntBoolPair(5,(matchHadp&&matchHadq)||(matchHadpq&&matchHadqp))); 

  //obs6 : b-tagging information
  double Obs6 = 0;
  if ( fabs(solution.getHadb().getBDiscriminator("trackCountingJetTags") +10) < 0.0001 || fabs(solution.getLepb().getBDiscriminator("trackCountingJetTags") +10)< 0.0001 ){
    Obs6 = -10.;
  } else {
    Obs6 = (solution.getHadb().getBDiscriminator("trackCountingJetTags")+solution.getLepb().getBDiscriminator("trackCountingJetTags"));
  }
  evtselectVarVal.push_back(IntDblPair(6,Obs6)); 
  evtselectVarMatch.push_back(IntBoolPair(6,1)); 
 
  //obs7 : chi2 value of kinematical fit with W-mass constraint
  double Obs7 =0;
  if(solution.getProbChi2() <0){Obs7 = -0;} else { Obs7 = log10(solution.getProbChi2()+.00001);}
  evtselectVarVal.push_back(IntDblPair(7,Obs7)); 
  evtselectVarMatch.push_back(IntBoolPair(7,((matchHadp&&matchHadq)||(matchHadpq&&matchHadqp)))); 

  //obs8(=7+1)
  double Obs8 =  solution.getCalHadt().p4().pt();
  evtselectVarVal.push_back(IntDblPair(8,Obs8));
  evtselectVarMatch.push_back(IntBoolPair(8, matchHadb&&((matchHadp&&matchHadq)||(matchHadpq&&matchHadqp))));

  //obs9
  double Obs9  = fabs(solution.getCalHadt().p4().eta());
  evtselectVarVal.push_back(IntDblPair(9,Obs9));
  evtselectVarMatch.push_back(IntBoolPair(9, matchHadb&&((matchHadp&&matchHadq)||(matchHadpq&&matchHadqp))));

  //obs10  
  double Obs10  = solution.getCalHadt().p4().theta();
  evtselectVarVal.push_back(IntDblPair(10,Obs10));
  evtselectVarMatch.push_back(IntBoolPair(10, matchHadb&&((matchHadp&&matchHadq)||(matchHadpq&&matchHadqp))));
  
  //obs11
  double Obs11  = solution.getCalHadW().p4().pt();
  evtselectVarVal.push_back(IntDblPair(11,Obs11));
  evtselectVarMatch.push_back(IntBoolPair(11, (matchHadp&&matchHadq)||(matchHadpq&&matchHadqp)));
  
  //obs12
  double Obs12  = fabs(solution.getCalHadW().p4().eta());
  evtselectVarVal.push_back(IntDblPair(12,Obs12));
  evtselectVarMatch.push_back(IntBoolPair(12, (matchHadp&&matchHadq)||(matchHadpq&&matchHadqp)));
  
  //obs13
  double Obs13  = solution.getCalHadW().p4().theta();
  evtselectVarVal.push_back(IntDblPair(13,Obs13));
  evtselectVarMatch.push_back(IntBoolPair(13, (matchHadp&&matchHadq)||(matchHadpq&&matchHadqp)));
	
  //obs14
  double Obs14  = solution.getCalHadb().p4().pt();
  evtselectVarVal.push_back(IntDblPair(14,Obs14));
  evtselectVarMatch.push_back(IntBoolPair(14, matchHadb));
	    
  //obs15
  double Obs15  = fabs(solution.getCalHadb().p4().eta());
  evtselectVarVal.push_back(IntDblPair(15,Obs15));
  evtselectVarMatch.push_back(IntBoolPair(15, matchHadb));
	    
  //obs16
  double Obs16  = solution.getCalHadb().p4().theta();
  evtselectVarVal.push_back(IntDblPair(16,Obs16));
  evtselectVarMatch.push_back(IntBoolPair(16, matchHadb));
	    
  //obs17
  double Obs17  = solution.getCalHadp().p4().pt();
  evtselectVarVal.push_back(IntDblPair(17,Obs17));
  evtselectVarMatch.push_back(IntBoolPair(17, matchHadp));
	    
  //obs18
  double Obs18  = fabs(solution.getCalHadp().p4().eta());
  evtselectVarVal.push_back(IntDblPair(18,Obs18));
  evtselectVarMatch.push_back(IntBoolPair(18, matchHadp));
	    
  //obs19
  double Obs19  = solution.getCalHadp().p4().theta();
  evtselectVarVal.push_back(IntDblPair(19,Obs19));
  evtselectVarMatch.push_back(IntBoolPair(19, matchHadp));
	    
  //obs20
  double Obs20  = solution.getCalHadq().p4().pt();
  evtselectVarVal.push_back(IntDblPair(20,Obs20));
  evtselectVarMatch.push_back(IntBoolPair(20, matchHadq));
	    
  //obs21
  double Obs21  = fabs(solution.getCalHadq().p4().eta());
  evtselectVarVal.push_back(IntDblPair(21,Obs21));
  evtselectVarMatch.push_back(IntBoolPair(21, matchHadq));
	    
  //obs22
  double Obs22  = solution.getCalHadq().p4().theta();
  evtselectVarVal.push_back(IntDblPair(22,Obs22));
  evtselectVarMatch.push_back(IntBoolPair(22, matchHadq));
	    
  //obs23
  double Obs23  = solution.getCalLept().p4().pt();
  evtselectVarVal.push_back(IntDblPair(23,Obs23));
  evtselectVarMatch.push_back(IntBoolPair(23, matchLepl&&matchLepn&&matchLepb));
	    
  //obs24
  double Obs24  = fabs(solution.getCalLept().p4().eta());
  evtselectVarVal.push_back(IntDblPair(24,Obs24));
  evtselectVarMatch.push_back(IntBoolPair(24, matchLepl&&matchLepn&&matchLepb));
	    
  //obs25
  double Obs25  = solution.getCalLept().p4().theta();
  evtselectVarVal.push_back(IntDblPair(25,Obs25));
  evtselectVarMatch.push_back(IntBoolPair(25, matchLepl&&matchLepn&&matchLepb));
	    
  //obs26
  double Obs26  = solution.getRecLepW().p4().pt(); 
  evtselectVarVal.push_back(IntDblPair(26,Obs26));
  evtselectVarMatch.push_back(IntBoolPair(26, matchLepl&&matchLepn));

  //obs27
  double Obs27  = fabs(solution.getRecLepW().p4().eta());  
  evtselectVarVal.push_back(IntDblPair(27,Obs27));
  evtselectVarMatch.push_back(IntBoolPair(27, matchLepl&&matchLepn));
	    
  //obs28
  double Obs28  = solution.getRecLepW().p4().theta(); 
  evtselectVarVal.push_back(IntDblPair(28,Obs28));
  evtselectVarMatch.push_back(IntBoolPair(28, matchLepl&&matchLepn));
  
  //obs29 
  double Obs29  = solution.getCalLepb().p4().pt();
  evtselectVarVal.push_back(IntDblPair(29,Obs29));
  evtselectVarMatch.push_back(IntBoolPair(29, matchLepb));
  
  //obs30
  double Obs30  = fabs(solution.getCalLepb().p4().eta());
  evtselectVarVal.push_back(IntDblPair(30,Obs30));
  evtselectVarMatch.push_back(IntBoolPair(30, matchLepb));
  
  //obs31
  double Obs31  = solution.getCalLepb().p4().theta();
  evtselectVarVal.push_back(IntDblPair(31,Obs31));
  evtselectVarMatch.push_back(IntBoolPair(31, matchLepb));
  
  //obs32
  double Obs32;   
  if (solution.getDecay() == "muon") Obs32 = solution.getRecLepm().p4().pt(); 
  if (solution.getDecay() == "electron") Obs32 = solution.getRecLepe().p4().pt(); 
  evtselectVarVal.push_back(IntDblPair(32,Obs32));
  evtselectVarMatch.push_back(IntBoolPair(32, matchLepl));
   
  //obs33
  double Obs33;
  if (solution.getDecay() == "muon") Obs33 = fabs(solution.getRecLepm().p4().eta());
  if (solution.getDecay() == "electron") Obs33 = fabs(solution.getRecLepe().p4().eta());
  evtselectVarVal.push_back(IntDblPair(33,Obs33));
  evtselectVarMatch.push_back(IntBoolPair(33, matchLepl));
 
  //obs34
  double Obs34;
  if (solution.getDecay() == "muon") Obs34 = fabs(solution.getRecLepm().p4().theta());
  if (solution.getDecay() == "electron") Obs34 = fabs(solution.getRecLepe().p4().theta()); 
  evtselectVarVal.push_back(IntDblPair(34,Obs34));
  evtselectVarMatch.push_back(IntBoolPair(34, matchLepl));
  
  //obs35
  double Obs35  = solution.getFitLepn().p4().pt();
  evtselectVarVal.push_back(IntDblPair(35,Obs35));
  evtselectVarMatch.push_back(IntBoolPair(35, matchLepn));
  
  //obs36
  double Obs36  = fabs(solution.getFitLepn().p4().eta());
  evtselectVarVal.push_back(IntDblPair(36,Obs36));
  evtselectVarMatch.push_back(IntBoolPair(36, matchLepn));
  
  //obs37
  double Obs37  = solution.getFitLepn().p4().theta();
  evtselectVarVal.push_back(IntDblPair(37,Obs37));
  evtselectVarMatch.push_back(IntBoolPair(37, matchLepn));
  
  // 2 particle kinematics
  //obs38 
  double Obs38  = fabs(solution.getCalHadW().p4().phi()- solution.getRecLepW().p4().phi());
  if (Obs38 > 3.1415927)  Obs38 =  2*3.1415927 - Obs31;
  if (Obs38 < -3.1415927) Obs38 = -2*3.1415927 - Obs31;
  evtselectVarVal.push_back(IntDblPair(38,Obs38));
  evtselectVarMatch.push_back(IntBoolPair(38, matchLepl&&matchLepn&&((matchHadp&&matchHadq)||(matchHadpq&&matchHadqp))));
  
  //obs39
  double Obs39  = fabs(solution.getCalHadW().p4().eta()- solution.getRecLepW().p4().eta());
  evtselectVarVal.push_back(IntDblPair(39,Obs39));
  evtselectVarMatch.push_back(IntBoolPair(39, matchLepl&&matchLepn&&((matchHadp&&matchHadq)||(matchHadpq&&matchHadqp))));
  
  //obs40
  double Obs40  = fabs(solution.getCalHadW().p4().theta()- solution.getRecLepW().p4().theta());
  evtselectVarVal.push_back(IntDblPair(40,Obs40));
  evtselectVarMatch.push_back(IntBoolPair(40, matchLepl&&matchLepn&&((matchHadp&&matchHadq)||(matchHadpq&&matchHadqp))));
  
  //obs41
  double Obs41  = ROOT::Math::VectorUtil::DeltaR(solution.getCalHadW().p4(), solution.getRecLepW().p4());
  evtselectVarVal.push_back(IntDblPair(41,Obs41));
  evtselectVarMatch.push_back(IntBoolPair(41, matchLepl&&matchLepn&&((matchHadp&&matchHadq)||(matchHadpq&&matchHadqp))));
  
  //obs42
  double Obs42  = fabs(solution.getCalHadb().p4().phi()- solution.getCalLepb().p4().phi());
  if (Obs42 > 3.1415927)  Obs42 =  2*3.1415927 - Obs42;
  if (Obs42 < -3.1415927) Obs42 = -2*3.1415927 - Obs42;
  evtselectVarVal.push_back(IntDblPair(42,Obs42));
  evtselectVarMatch.push_back(IntBoolPair(42, matchHadb&&matchLepb));
  
  //obs43
  double Obs43  = fabs(solution.getCalHadb().p4().eta()- solution.getCalLepb().p4().eta());
  evtselectVarVal.push_back(IntDblPair(43,Obs43));
  evtselectVarMatch.push_back(IntBoolPair(43, matchHadb&&matchLepb));
  
  //obs44
  double Obs44 = fabs(solution.getCalHadb().p4().theta()- solution.getCalLepb().p4().theta());
  evtselectVarVal.push_back(IntDblPair(44,Obs44));
  evtselectVarMatch.push_back(IntBoolPair(44, matchHadb&&matchLepb));
  
  //obs45
  double Obs45  = ROOT::Math::VectorUtil::DeltaR(solution.getCalHadb().p4(), solution.getCalLepb().p4());
  evtselectVarVal.push_back(IntDblPair(45,Obs45));
  evtselectVarMatch.push_back(IntBoolPair(45, matchHadb&&matchLepb));
  
  //obs46
  double Obs46  = fabs(solution.getCalHadb().p4().phi()- solution.getCalHadW().p4().phi());
  if (Obs46 > 3.1415927)  Obs46 =  2*3.1415927 - Obs46;
  if (Obs46 < -3.1415927) Obs46 = -2*3.1415927 - Obs46;
  evtselectVarVal.push_back(IntDblPair(46,Obs46));
  evtselectVarMatch.push_back(IntBoolPair(46, matchHadb&&((matchHadq&&matchHadp)||(matchHadpq&&matchHadqp))));

  //obs47
  double Obs47  = fabs(solution.getCalHadb().p4().eta()- solution.getCalHadW().p4().eta());
  evtselectVarVal.push_back(IntDblPair(47,Obs47));
  evtselectVarMatch.push_back(IntBoolPair(47, matchHadb&&((matchHadp&&matchHadq)||(matchHadpq&&matchHadqp))));
  
  //obs48
  double Obs48  = fabs(solution.getCalHadb().p4().theta()- solution.getCalHadW().p4().theta());
  evtselectVarVal.push_back(IntDblPair(48,Obs48));
  evtselectVarMatch.push_back(IntBoolPair(48, matchHadb&&((matchHadp&&matchHadq)||(matchHadpq&&matchHadqp))));
  
  //obs49
  double Obs49  = ROOT::Math::VectorUtil::DeltaR(solution.getCalHadb().p4(), solution.getCalHadW().p4());
  evtselectVarVal.push_back(IntDblPair(49,Obs49));
  evtselectVarMatch.push_back(IntBoolPair(49, matchHadb&&((matchHadp&&matchHadq)||(matchHadpq&&matchHadqp))));
  
  //obs50
  double Obs50  = fabs(solution.getCalLepb().p4().phi()- solution.getRecLepW().p4().phi());
  if (Obs50 > 3.1415927)  Obs50 =  2*3.1415927 - Obs50;
  if (Obs50 < -3.1415927) Obs50 = -2*3.1415927 - Obs50;
  evtselectVarVal.push_back(IntDblPair(50,Obs50));
  evtselectVarMatch.push_back(IntBoolPair(50, matchLepb&&matchLepn&&matchLepl));
  
  //obs51
  double Obs51  = fabs(solution.getCalLepb().p4().eta()- solution.getRecLepW().p4().eta()); 
  evtselectVarVal.push_back(IntDblPair(51,Obs51));
  evtselectVarMatch.push_back(IntBoolPair(51, matchLepb&&matchLepn&&matchLepl));
  
  //obs52
  double Obs52  = fabs(solution.getCalLepb().p4().theta()- solution.getRecLepW().p4().theta());
  evtselectVarVal.push_back(IntDblPair(52,Obs52));
  evtselectVarMatch.push_back(IntBoolPair(52, matchLepb&&matchLepn&&matchLepl));
  
  //obs53
  double Obs53  = ROOT::Math::VectorUtil::DeltaR(solution.getCalLepb().p4(), solution.getRecLepW().p4());
  evtselectVarVal.push_back(IntDblPair(53,Obs53));
  evtselectVarMatch.push_back(IntBoolPair(53, matchLepb&&matchLepn&&matchLepl));
  
  //obs54
  double Obs54 = fabs(solution.getCalHadp().p4().phi()- solution.getCalHadq().p4().phi());
  if (Obs54 > 3.1415927)  Obs54 =  2*3.1415927 - Obs54;
  if (Obs54 < -3.1415927) Obs54 = -2*3.1415927 - Obs54;
  evtselectVarVal.push_back(IntDblPair(54,Obs54));
  evtselectVarMatch.push_back(IntBoolPair(54, (matchHadp&&matchHadq)||(matchHadpq&&matchHadqp)));

  //obs55
  double Obs55 = fabs(solution.getCalHadp().p4().eta()- solution.getCalHadq().p4().eta());
  evtselectVarVal.push_back(IntDblPair(55,Obs55));
  evtselectVarMatch.push_back(IntBoolPair(55, (matchHadp&&matchHadq)||(matchHadpq&&matchHadqp)));
  
  //obs56
  double Obs56  = fabs(solution.getCalHadp().p4().theta()- solution.getCalHadq().p4().theta());
  evtselectVarVal.push_back(IntDblPair(56,Obs56));
  evtselectVarMatch.push_back(IntBoolPair(56, (matchHadp&&matchHadq)||(matchHadpq&&matchHadqp)));
  
  //obs57
  double Obs57  = ROOT::Math::VectorUtil::DeltaR(solution.getCalHadp().p4(), solution.getCalHadq().p4());
  evtselectVarVal.push_back(IntDblPair(57,Obs57));
  evtselectVarMatch.push_back(IntBoolPair(57, (matchHadp&&matchHadq)||(matchHadpq&&matchHadqp)));
  
  //obs58
  double Obs58;
  if (solution.getDecay() == "muon") Obs58  = fabs(solution.getRecLepm().p4().phi()- solution.getRecLepn().p4().phi()); 
  if (solution.getDecay() == "electron") Obs58 = fabs(solution.getRecLepe().p4().phi()- solution.getRecLepn().p4().phi()); 
  if (Obs58 > 3.1415927)  Obs58 =  2*3.1415927 - Obs58;
  if (Obs58 < -3.1415927) Obs58 = -2*3.1415927 - Obs58;
  evtselectVarVal.push_back(IntDblPair(58,Obs58));
  evtselectVarMatch.push_back(IntBoolPair(58, matchLepl&&matchLepn));
 
  //obs59
  double Obs59;
  if (solution.getDecay() == "muon") Obs59 = fabs(solution.getRecLepm().p4().eta()- solution.getRecLepn().p4().eta()); 
  if (solution.getDecay() == "electron") Obs59 = fabs(solution.getRecLepe().p4().eta()- solution.getRecLepn().p4().eta());  
  evtselectVarVal.push_back(IntDblPair(59,Obs59));
  evtselectVarMatch.push_back(IntBoolPair(59, matchLepl&&matchLepn));
  
  //obs60
  double Obs60;
  if (solution.getDecay() == "muon") Obs60  = fabs(solution.getRecLepm().p4().theta()- solution.getRecLepn().p4().theta());
  if (solution.getDecay() == "electron") Obs60  = fabs(solution.getRecLepe().p4().theta()- solution.getRecLepn().p4().theta());
  evtselectVarVal.push_back(IntDblPair(60,Obs60));
  evtselectVarMatch.push_back(IntBoolPair(60, matchLepl&&matchLepn));
  
  //obs61 
  double Obs61; 
  if (solution.getDecay() == "muon") Obs61 = ROOT::Math::VectorUtil::DeltaR(solution.getRecLepm().p4(), solution.getRecLepn().p4());
  if (solution.getDecay() == "electron") Obs61 = ROOT::Math::VectorUtil::DeltaR(solution.getRecLepe().p4(), solution.getRecLepn().p4());
  evtselectVarVal.push_back(IntDblPair(61,Obs61));
  evtselectVarMatch.push_back(IntBoolPair(61, matchLepl&&matchLepn));
  
  // miscellaneous event variables
   

  //obs62
  double Obs62  = (jets->size() > 4 ? (*jets)[4].p4().Et() / (*jets)[3].p4().Et() : 0);
  //double Obs62 = 1;
  evtselectVarVal.push_back(IntDblPair(62,Obs62));
  evtselectVarMatch.push_back(IntBoolPair(62, ((matchHadp&&matchHadq)||(matchHadpq&&matchHadqp))&&matchHadb&&matchLepb));
  
  float calJetsSumEt = 0;
  for (unsigned int i = 4; i < jets->size(); i++) {
    calJetsSumEt += (*jets)[i].p4().Et();
  }
  
  //obs63
  double Obs63  = calJetsSumEt / ((*jets)[0].p4().Et()+(*jets)[1].p4().Et()+(*jets)[2].p4().Et()+(*jets)[3].p4().Et()+(*jets)[4].p4().Et());
  //double Obs63 =1;  
  evtselectVarVal.push_back(IntDblPair(63,Obs63));
  evtselectVarMatch.push_back(IntBoolPair(63, ((matchHadp&&matchHadq)||(matchHadpq&&matchHadqp))&&matchHadb&&matchLepb));
  
  //obs64
  double Obs64  = solution.getProbChi2();
  evtselectVarVal.push_back(IntDblPair(64,Obs64));
  evtselectVarMatch.push_back(IntBoolPair(64, ((matchHadp&&matchHadq)||(matchHadpq&&matchHadqp))&&matchHadb&&matchLepb));
  
  //obs65
  double Obs65  = solution.getFitHadt().p4().mass() - solution.getCalHadt().p4().mass();
  evtselectVarVal.push_back(IntDblPair(65,Obs65));
  evtselectVarMatch.push_back(IntBoolPair(65, ((matchHadp&&matchHadq)||(matchHadpq&&matchHadqp))&&matchHadb&&matchLepb)); 
  
  //obs66
  double Obs66  = solution.getFitLept().p4().mass() - solution.getCalLept().p4().mass();
  evtselectVarVal.push_back(IntDblPair(66,Obs66));
  evtselectVarMatch.push_back(IntBoolPair(66, ((matchHadp&&matchHadq)||(matchHadpq&&matchHadqp))&&matchHadb&&matchLepb)); 
  
  if(debug) cout << "observables calculated" << endl;

  if (!matchOnly) solution.setLRJetCombObservables(evtselectVarVal);
  return evtselectVarMatch;
}

#include "TopQuarkAnalysis/TopJetCombination/plugins/TtFullLepHypGenMatch.h"
#include "AnalysisDataFormats/TopObjects/interface/TtFullLepEvtPartons.h"
#include "DataFormats/Math/interface/deltaR.h"

TtFullLepHypGenMatch::TtFullLepHypGenMatch(const edm::ParameterSet& cfg):
  TtFullLepHypothesis( cfg ) 
{  
}

TtFullLepHypGenMatch::~TtFullLepHypGenMatch() { }

void
TtFullLepHypGenMatch::buildHypo(edm::Event& evt,
			        const edm::Handle<std::vector<pat::Electron > >& elecs, 
			        const edm::Handle<std::vector<pat::Muon> >& mus, 
			        const edm::Handle<std::vector<pat::Jet> >& jets, 
			        const edm::Handle<std::vector<pat::MET> >& mets, 
			        std::vector<int>& match,
				const unsigned int iComb)
{
  // -----------------------------------------------------
  // add jets
  // -----------------------------------------------------
  for(unsigned idx=0; idx<match.size(); ++idx){
    if( isValid(match[idx], jets) ){
      switch(idx){
      case TtFullLepEvtPartons::B:
	setCandidate(jets, match[idx], b_); break;
      case TtFullLepEvtPartons::BBar:
	setCandidate(jets, match[idx], bBar_); break;	
      }
    }
  }

  // -----------------------------------------------------
  // add leptons
  // -----------------------------------------------------
  if( !mus->empty() ){
    std::pair<int, int> ijLeptons = findMatchingLeptons(evt,mus);
    if( ijLeptons.first >=0 )
      setCandidate(mus, ijLeptons.first, lepton_);
    match.push_back( ijLeptons.first );
    if( ijLeptons.second >=0 )
      setCandidate(mus, ijLeptons.second, leptonBar_);
    match.push_back( ijLeptons.second );    
  }
  else{
    match.push_back( -1 );
    match.push_back( -1 );
  }
  // -----------------------------------------------------
  // add met and neutrinos
  // -----------------------------------------------------  
  if( !mets->empty() ){
    //setCandidate(mets, 0, met_);
    buildMatchingNeutrinos(evt, mets);  
  }    
}

std::pair<int, int>
TtFullLepHypGenMatch::findMatchingLeptons(edm::Event& evt, const edm::Handle<std::vector<pat::Muon> >& leps)
{
  std::pair<int, int> genIdcs(-1,-1);

  // get genEvent
  edm::Handle<TtGenEvent> genEvt;
  evt.getByLabel("genEvt", genEvt);  
  
  if( genEvt->isTtBar() && genEvt->isFullLeptonic() && genEvt->lepton() && genEvt->leptonBar() ){
    double minDRlep = -1;
    double minDRlepbar = -1;
    for(unsigned i=0; i<leps->size(); ++i){
      double dRlep = deltaR(genEvt->lepton()->eta(), genEvt->lepton()->phi(), (*leps)[i].eta(), (*leps)[i].phi());
      if(minDRlep<0 || dRlep<minDRlep){
	minDRlep=dRlep;
	genIdcs.first=i;
      }
      double dRlepbar = deltaR(genEvt->leptonBar()->eta(), genEvt->leptonBar()->phi(), (*leps)[i].eta(), (*leps)[i].phi());
      if(minDRlepbar<0 || dRlepbar<minDRlepbar){
	minDRlepbar=dRlepbar;
	genIdcs.second=i;
      }      
    }
  }
  return genIdcs;
}

void
TtFullLepHypGenMatch::buildMatchingNeutrinos(edm::Event& evt, const edm::Handle<std::vector<pat::MET> >& mets)
{
  // get genEvent
  edm::Handle<TtGenEvent> genEvt;
  evt.getByLabel("genEvt", genEvt);
  
  if( genEvt->isTtBar() && genEvt->isFullLeptonic() && genEvt->neutrino() && genEvt->neutrinoBar() ){
    double momXNu    = genEvt->neutrino()   ->px();
    double momYNu    = genEvt->neutrino()   ->py(); 
    double momXNuBar = genEvt->neutrinoBar()->px();
    double momYNuBar = genEvt->neutrinoBar()->py();
        
    double momXMet = mets->at(0).px();
    double momYMet = mets->at(0).py();

    double momXNeutrino = 0.5*(momXNu - momXNuBar + momXMet);
    double momYNeutrino = 0.5*(momYNu - momYNuBar + momYMet);   
    double momXNeutrinoBar = momXMet - momXNeutrino;
    double momYNeutrinoBar = momYMet - momYNeutrino; 
  
    math::XYZTLorentzVector recNuFM(momXNeutrino,momYNeutrino,0,sqrt(momXNeutrino * momXNeutrino + momYNeutrino * momYNeutrino));    
    recNu = new reco::LeafCandidate(0,recNuFM);
    
    math::XYZTLorentzVector recNuBarFM(momXNeutrinoBar,momYNeutrinoBar,0,sqrt(momXNeutrinoBar * momXNeutrinoBar + momYNeutrinoBar * momYNeutrinoBar));
    recNuBar = new reco::LeafCandidate(0,recNuBarFM);		      
  }
}

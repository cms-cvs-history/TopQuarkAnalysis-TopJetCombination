#include "TopQuarkAnalysis/TopJetCombination/plugins/TtSemiHypothesisGenMatch.h"
#include "TopQuarkAnalysis/TopTools/interface/TtSemiEvtPartons.h"
#include "DataFormats/Math/interface/deltaR.h"

TtSemiHypothesisGenMatch::TtSemiHypothesisGenMatch(const edm::ParameterSet& cfg):
  TtSemiHypothesis( cfg ) { }

TtSemiHypothesisGenMatch::~TtSemiHypothesisGenMatch() { }

void
TtSemiHypothesisGenMatch::buildHypo(edm::Event& evt,
				    const edm::Handle<edm::View<reco::RecoCandidate> >& leps, 
				    const edm::Handle<std::vector<pat::MET> >& mets, 
				    const edm::Handle<std::vector<pat::Jet> >& jets, 
				    std::vector<int>& match)
{
  // -----------------------------------------------------
  // add jets
  // -----------------------------------------------------
  for(unsigned idx=0; idx<match.size(); ++idx){
    if( isValid(match[idx], jets) ){
      switch(idx){
      case TtSemiEvtPartons::LightQ:
	setCandidate(jets, match[idx], lightQ_); break;
      case TtSemiEvtPartons::LightQBar:
	setCandidate(jets, match[idx], lightQBar_); break;
      case TtSemiEvtPartons::HadB:
	setCandidate(jets, match[idx], hadronicB_); break;
      case TtSemiEvtPartons::LepB: 
	setCandidate(jets, match[idx], leptonicB_); break;
      }
    }
  }

  // -----------------------------------------------------
  // add lepton
  // -----------------------------------------------------
  if( !leps->empty() ){
    int iLepton = findMatchingLepton(evt,leps);
    if( iLepton>=0 )
      setCandidate(leps, iLepton, lepton_);
  }
  
  // -----------------------------------------------------
  // add neutrino
  // -----------------------------------------------------
  if( !mets->empty() )
    setCandidate(mets, 0, neutrino_);
}

int
TtSemiHypothesisGenMatch::findMatchingLepton(edm::Event& evt, const edm::Handle<edm::View<reco::RecoCandidate> >& leps)
{
  int genIdx=-1;

  // set genEvent
  edm::Handle<TtGenEvent> genEvt;
  evt.getByLabel("genEvt", genEvt);  
  
  if( genEvt->isTtBar() && genEvt->isSemiLeptonic() && genEvt->lepton() ){
    double minDR=-1;
    for(unsigned i=0; i<leps->size(); ++i){
      double dR = deltaR(genEvt->lepton()->eta(), genEvt->lepton()->phi(), (*leps)[i].eta(), (*leps)[i].phi());
      if(minDR<0 || dR<minDR){
	minDR=dR;
	genIdx=i;
      }
    }
  }
  return genIdx;
}

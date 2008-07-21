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
  // add jets; the order of match is Q, QBar, hadB, lepB
  // -----------------------------------------------------
  for(unsigned idx=0; idx<match.size(); ++idx){    
    if( isValid(match[idx], jets) ){
      edm::Ref<std::vector<pat::Jet> > ref=edm::Ref<std::vector<pat::Jet> >(jets, match[idx]);
      reco::ShallowCloneCandidate buffer(reco::CandidateBaseRef( ref ), ref->charge(), ref->p4(), ref->vertex());
      switch(idx){
      case TtSemiEvtPartons::LightQ: 
	lightQ_   = new reco::ShallowCloneCandidate( buffer ); break;
      case TtSemiEvtPartons::LightQBar: 
	lightQBar_= new reco::ShallowCloneCandidate( buffer ); break;
      case TtSemiEvtPartons::HadB: 
	hadronicB_= new reco::ShallowCloneCandidate( buffer ); break;
      case TtSemiEvtPartons::LepB: 
	leptonicB_= new reco::ShallowCloneCandidate( buffer ); break;
      }
    }
  }

  // -----------------------------------------------------
  // add lepton
  // -----------------------------------------------------
  if( !leps->empty() ){
    int iLepton = findMatchingLepton(evt,leps);
    if( iLepton>=0 ){
      edm::Ref<edm::View<reco::RecoCandidate> > ref=edm::Ref<edm::View<reco::RecoCandidate> >(leps, iLepton);
      reco::ShallowCloneCandidate buffer(reco::CandidateBaseRef( ref ), ref->charge(), ref->p4(), ref->vertex());
      lepton_= new reco::ShallowCloneCandidate( buffer );
    }
  }
  
  // -----------------------------------------------------
  // add neutrino
  // -----------------------------------------------------
  {
    if( !mets->empty() ){
      edm::Ref<std::vector<pat::MET> > ref=edm::Ref<std::vector<pat::MET> >(mets, 0);
      reco::ShallowCloneCandidate buffer(reco::CandidateBaseRef( ref ), ref->charge(), ref->p4(), ref->vertex());
      neutrino_= new reco::ShallowCloneCandidate( buffer );
    }
  }
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

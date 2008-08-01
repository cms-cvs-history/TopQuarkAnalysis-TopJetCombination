#include "TopQuarkAnalysis/TopJetCombination/interface/TtSemiHypothesis.h"

#include "PhysicsTools/CandUtils/interface/AddFourMomenta.h"

TtSemiHypothesis::TtSemiHypothesis(const edm::ParameterSet& cfg):
  jets_(cfg.getParameter<edm::InputTag>("jets")),
  leps_(cfg.getParameter<edm::InputTag>("leps")),
  mets_(cfg.getParameter<edm::InputTag>("mets")),
  lightQ_(0), lightQBar_(0), hadronicB_(0), 
  leptonicB_(0), neutrino_(0), lepton_(0)
{
  getMatch_ = false;
  if(cfg.exists("match")) {
    getMatch_ = true;
    match_ = cfg.getParameter<edm::InputTag>("match");
  }

  produces<reco::NamedCompositeCandidate>();
  produces<int>("Key");
  produces<std::vector<int> >("Match");
}

TtSemiHypothesis::~TtSemiHypothesis()
{
  if( lightQ_   ) delete lightQ_;
  if( lightQBar_) delete lightQBar_;
  if( hadronicB_) delete hadronicB_;
  if( leptonicB_) delete leptonicB_;
  if( neutrino_ ) delete neutrino_;
  if( lepton_   ) delete lepton_;
}

void
TtSemiHypothesis::produce(edm::Event& evt, const edm::EventSetup& setup)
{
  edm::Handle<std::vector<pat::Jet> > jets;
  evt.getByLabel(jets_, jets);
  
  edm::Handle<edm::View<reco::RecoCandidate> > leps;
  evt.getByLabel(leps_, leps);

  edm::Handle<std::vector<pat::MET> > mets;
  evt.getByLabel(mets_, mets);

  std::vector<int> match;
  if(getMatch_) {
    edm::Handle<std::vector<int> > matchHandle;
    evt.getByLabel(match_, matchHandle);
    match = *matchHandle;
  }

  // feed out hyp
  std::auto_ptr<reco::NamedCompositeCandidate> pOut(new reco::NamedCompositeCandidate);
  buildHypo(evt, leps, mets, jets, match);
  *pOut=hypo();
  evt.put(pOut);
  
  // feed out key
  std::auto_ptr<int> pKey(new int);
  buildKey();
  *pKey=key();
  evt.put(pKey, "Key");

  // feed out match
  std::auto_ptr<std::vector<int> > pMatch(new std::vector<int>);
  for(unsigned int i=0; i<match.size(); ++i)
    pMatch->push_back( match[i] );
  evt.put(pMatch, "Match");
}

void
TtSemiHypothesis::setCandidate(const edm::Handle<pat::Particle>& handle, reco::ShallowCloneCandidate* &clone)
{
  edm::RefProd<pat::Particle> ref=edm::RefProd<pat::Particle>( handle );
  reco::ShallowCloneCandidate buffer(reco::CandidateBaseRef( ref ), ref->charge(), ref->p4(), ref->vertex());
  clone = new reco::ShallowCloneCandidate( buffer );  
}

reco::NamedCompositeCandidate
TtSemiHypothesis::hypo()
{
  // -------------------------------------------------------------------------
  // note that there is a persistency issue with 1612 were reflex classes are 
  // inadequate for reco::RecoCandidate which is in the edm::View the generic 
  // lepton. TtSemiHypotheses may only be made persistent when commenting the 
  // leptonic branch o0f the hypothesis. This problem does not exist any more 
  // in the 21X series
  // -------------------------------------------------------------------------

  // check for sanity of the hypothesis
  if( !lightQ_ || !lightQBar_ || !hadronicB_ || 
      !leptonicB_ ||  !neutrino_ || !lepton_ )
    return reco::NamedCompositeCandidate();
  
  // setup transient references
  reco::NamedCompositeCandidate hyp, hadTop, hadW, lepTop, lepW;

  AddFourMomenta addFourMomenta;  
  // build up the top branch that decays leptonically
  lepW  .addDaughter(*lepton_,   TtSemiDaughter::Lep    );
  lepW  .addDaughter(*neutrino_, TtSemiDaughter::Nu     );
  addFourMomenta.set( lepW );
  lepTop.addDaughter( lepW,      TtSemiDaughter::LepW   );
  lepTop.addDaughter(*leptonicB_,TtSemiDaughter::LepB   );
  addFourMomenta.set( lepTop );
  
  // build up the top branch that decays hadronically
  hadW  .addDaughter(*lightQ_,   TtSemiDaughter::HadQ   );
  hadW  .addDaughter(*lightQBar_,TtSemiDaughter::HadP   );
  addFourMomenta.set( hadW );
  hadTop.addDaughter( hadW,      TtSemiDaughter::HadW   );
  hadTop.addDaughter(*hadronicB_,TtSemiDaughter::HadB   );
  addFourMomenta.set( hadTop );
  
  // build ttbar hypotheses
  hyp.addDaughter( lepTop,       TtSemiDaughter::LepTop );
  hyp.addDaughter( hadTop,       TtSemiDaughter::HadTop );
  addFourMomenta.set( hyp );
  
  return hyp;
}

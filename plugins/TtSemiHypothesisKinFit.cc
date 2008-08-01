#include "TopQuarkAnalysis/TopJetCombination/plugins/TtSemiHypothesisKinFit.h"

TtSemiHypothesisKinFit::TtSemiHypothesisKinFit(const edm::ParameterSet& cfg):
  TtSemiHypothesis( cfg ),
  status_(cfg.getParameter<edm::InputTag>("status")),
  hadB_  (cfg.getParameter<edm::InputTag>("hadB")),
  hadP_  (cfg.getParameter<edm::InputTag>("hadP")),
  hadQ_  (cfg.getParameter<edm::InputTag>("hadQ")),
  lepB_  (cfg.getParameter<edm::InputTag>("lepB")),
  lepL_  (cfg.getParameter<edm::InputTag>("lepL")),
  lepN_  (cfg.getParameter<edm::InputTag>("lepN"))
{
}

TtSemiHypothesisKinFit::~TtSemiHypothesisKinFit() { }

void
TtSemiHypothesisKinFit::buildHypo(edm::Event& evt,
				  const edm::Handle<edm::View<reco::RecoCandidate> >& leps, 
				  const edm::Handle<std::vector<pat::MET> >& mets, 
				  const edm::Handle<std::vector<pat::Jet> >& jets, 
				  std::vector<int>& match)
{
  edm::Handle<int> status;
  evt.getByLabel(status_, status);
  if(*status!=0){
    // create empty hypothesis if kinematic fit did not converge
    return;
  }

  // -----------------------------------------------------
  // add jets
  // -----------------------------------------------------
  edm::Handle<pat::Particle> hadB;
  edm::Handle<pat::Particle> hadP;
  edm::Handle<pat::Particle> hadQ;
  edm::Handle<pat::Particle> lepB;
  evt.getByLabel(hadB_, hadB);
  evt.getByLabel(hadP_, hadP);
  evt.getByLabel(hadQ_, hadQ);
  evt.getByLabel(lepB_, lepB);
  setCandidate(hadB, hadronicB_);
  setCandidate(hadP, lightQ_   );
  setCandidate(hadQ, lightQBar_);
  setCandidate(lepB, leptonicB_);

  // -----------------------------------------------------
  // add lepton
  // -----------------------------------------------------
  edm::Handle<pat::Particle> lepL;
  evt.getByLabel(lepL_, lepL);
  setCandidate(lepL, lepton_);

  // -----------------------------------------------------
  // add neutrino
  // -----------------------------------------------------
  edm::Handle<pat::Particle> lepN;
  evt.getByLabel(lepN_, lepN);
  setCandidate(lepN, neutrino_);
}

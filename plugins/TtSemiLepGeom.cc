#include "TopQuarkAnalysis/TopJetCombination/plugins/TtSemiLepGeom.h"
#include "TopQuarkAnalysis/TopTools/interface/TtSemiEvtPartons.h"

#include <Math/VectorUtil.h>

TtSemiLepGeom::TtSemiLepGeom(const edm::ParameterSet& cfg):
  TtSemiLepHypothesis( cfg ),  
  maxNJets_ (cfg.getParameter<int> ("maxNJets" )),
  useDeltaR_(cfg.getParameter<bool>("useDeltaR"))
{ }

TtSemiLepGeom::~TtSemiLepGeom() { }

void
TtSemiLepGeom::buildHypo(edm::Event& evt,
				const edm::Handle<edm::View<reco::RecoCandidate> >& leps,
				const edm::Handle<std::vector<pat::MET> >& mets, 
				const edm::Handle<std::vector<pat::Jet> >& jets, 
				std::vector<int>& match)
{
  if(leps->empty() || mets->empty() || (int)(jets->size())<maxNJets_ || maxNJets_<4){
    // create empty hypothesis
    return;
  }

  for(unsigned int i=0; i<4; ++i)
    match.push_back(-1);
  
  // -----------------------------------------------------
  // associate those two jets to the hadronic W boson that
  // have the smallest distance to each other
  // -----------------------------------------------------
  double minDist=-1.;
  int lightQ   =-1;
  int lightQBar=-1;
  for(int idx=0; idx<maxNJets_; ++idx){
    for(int jdx=(idx+1); jdx<maxNJets_; ++jdx){
      double dist = distance((*jets)[idx].p4(), (*jets)[jdx].p4());
      if( minDist<0. || dist<minDist ){
	minDist=dist;
	lightQ   =idx;
	lightQBar=jdx;
      }
    }
  }

  reco::Particle::LorentzVector wHad = (*jets)[lightQ].p4() + (*jets)[lightQBar].p4();

  // -----------------------------------------------------
  // associate to the hadronic b quark the remaining jet
  // that has the smallest distance to the hadronic W 
  // -----------------------------------------------------
  minDist=-1.;
  int hadB=-1;
  for(int idx=0; idx<maxNJets_; ++idx){
    // make sure it's not used up already from the hadronic W
    if( idx!=lightQ && idx!=lightQBar ) {
      double dist = distance((*jets)[idx].p4(), wHad);
      if( minDist<0. || dist<minDist ){
	minDist=dist;
	hadB=idx;
      }
    }
  }

  // -----------------------------------------------------
  // associate to the leptonic b quark the remaining jet
  // that has the smallest distance to the leading lepton
  // -----------------------------------------------------
  minDist=-1.;
  int lepB=-1;
  for(int idx=0; idx<maxNJets_; ++idx){
    // make sure it's not used up already from the hadronic decay chain
    if( idx!=lightQ && idx!=lightQBar && idx!=hadB ){
      double dist = distance((*jets)[idx].p4(), (*leps)[0].p4());
      if( minDist<0. || dist<minDist ){
	minDist=dist;
	lepB=idx;
      }
    }
  }
  
  // -----------------------------------------------------
  // add jets
  // -----------------------------------------------------
  if( isValid(lightQ, jets) ){
    setCandidate(jets, lightQ, lightQ_);
    match[TtSemiEvtPartons::LightQ] = lightQ;
  }

  if( isValid(lightQBar, jets) ){
    setCandidate(jets, lightQBar, lightQBar_);
    match[TtSemiEvtPartons::LightQBar] = lightQBar;
  }

  if( isValid(hadB, jets) ){
    setCandidate(jets, hadB, hadronicB_);
    match[TtSemiEvtPartons::HadB] = hadB;
  }
  
  if( isValid(lepB, jets) ){
    setCandidate(jets, lepB, leptonicB_);
    match[TtSemiEvtPartons::LepB] = lepB;
  }

  // -----------------------------------------------------
  // add lepton
  // -----------------------------------------------------
  setCandidate(leps, 0, lepton_);
  
  // -----------------------------------------------------
  // add neutrino
  // -----------------------------------------------------
  setCandidate(mets, 0, neutrino_);
}

double
TtSemiLepGeom::distance(const math::XYZTLorentzVector& v1, const math::XYZTLorentzVector& v2)
{
  // calculate the distance between two lorentz vectors 
  // using DeltaR or DeltaTheta
  if(useDeltaR_) return ROOT::Math::VectorUtil::DeltaR(v1, v2);
  return fabs(v1.theta() - v2.theta());
}

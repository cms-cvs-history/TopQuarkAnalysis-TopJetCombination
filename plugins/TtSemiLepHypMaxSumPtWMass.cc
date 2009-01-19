#include "TopQuarkAnalysis/TopJetCombination/plugins/TtSemiLepHypMaxSumPtWMass.h"
#include "TopQuarkAnalysis/TopTools/interface/TtSemiLepEvtPartons.h"

TtSemiLepHypMaxSumPtWMass::TtSemiLepHypMaxSumPtWMass(const edm::ParameterSet& cfg):
  TtSemiLepHypothesis( cfg ),
  maxNJets_(cfg.getParameter<int>   ("maxNJets")),
  wMass_   (cfg.getParameter<double>("wMass"   ))
{
  if(maxNJets_<4 && maxNJets_!=-1)
    throw cms::Exception("WrongConfig") 
      << "Parameter maxNJets can not be set to " << maxNJets_ << ". \n"
      << "It has to be larger than 4 or can be set to -1 to take all jets.";
}

TtSemiLepHypMaxSumPtWMass::~TtSemiLepHypMaxSumPtWMass() { }

void
TtSemiLepHypMaxSumPtWMass::buildHypo(edm::Event& evt,
				     const edm::Handle<edm::View<reco::RecoCandidate> >& leps, 
				     const edm::Handle<std::vector<pat::MET> >& mets, 
				     const edm::Handle<std::vector<pat::Jet> >& jets, 
				     std::vector<int>& match, const unsigned int iComb)
{
  if(leps->empty() || mets->empty() || jets->size()<4){
    // create empty hypothesis
    return;
  }

  unsigned maxNJets = maxNJets_;
  if(maxNJets_ == -1 || (int)jets->size() < maxNJets_) maxNJets = jets->size();

  match.clear();
  for(unsigned int i=0; i<4; ++i)
    match.push_back(-1);

  // -----------------------------------------------------
  // associate those jets with maximum pt of the vectorial 
  // sum to the hadronic decay chain
  // -----------------------------------------------------
  double maxPt=-1.;
  std::vector<unsigned> maxPtIndices;
  for(unsigned idx=0; idx<maxNJets; ++idx){
    for(unsigned jdx=(idx+1); jdx<maxNJets; ++jdx){
      for(unsigned kdx=(jdx+1); kdx<maxNJets; ++kdx){
	reco::Particle::LorentzVector sum = 
	  (*jets)[idx].p4()+
	  (*jets)[jdx].p4()+
	  (*jets)[kdx].p4();
	if( maxPt<0. || maxPt<sum.pt() ){
	  maxPt=sum.pt();
	  maxPtIndices.clear();
	  maxPtIndices.push_back(idx);
	  maxPtIndices.push_back(jdx);
	  maxPtIndices.push_back(kdx);
	}
      }
    }
  }

  // -----------------------------------------------------
  // associate those jets that get closest to the W mass
  // with their invariant mass to the W boson
  // -----------------------------------------------------
  double wDist =-1.;
  std::vector<unsigned> closestToWMassIndices;
  for(unsigned idx=0; idx<maxPtIndices.size(); ++idx){  
    for(unsigned jdx=0; jdx<maxPtIndices.size(); ++jdx){  
      if( jdx==idx || maxPtIndices[idx]>maxPtIndices[jdx] ) continue;
      reco::Particle::LorentzVector sum = 
	(*jets)[maxPtIndices[idx]].p4()+
	(*jets)[maxPtIndices[jdx]].p4();
      if( wDist<0. || wDist>fabs(sum.mass()-wMass_) ){
	wDist=fabs(sum.mass()-wMass_);
	closestToWMassIndices.clear();
	closestToWMassIndices.push_back(maxPtIndices[idx]);
	closestToWMassIndices.push_back(maxPtIndices[jdx]);
      }
    }
  }

  // -----------------------------------------------------
  // associate the remaining jet with maximum pt of the   
  // vectorial sum with the leading lepton with the 
  // leptonic decay chain
  // -----------------------------------------------------
  maxPt=-1.;
  unsigned lepB=0;
  for(unsigned idx=0; idx<maxNJets; ++idx){
    // make sure it's not used up already from the hadronic decay chain
    if( std::find(maxPtIndices.begin(), maxPtIndices.end(), idx) == maxPtIndices.end() ){
      reco::Particle::LorentzVector sum = 
	(*jets)[idx].p4()+(*leps)[ 0 ].p4();
      if( maxPt<0. || maxPt<sum.pt() ){
	maxPt=sum.pt();
	lepB=idx;
      }
    }
  }

  // -----------------------------------------------------
  // add jets
  // -----------------------------------------------------
  if( isValid(closestToWMassIndices[0], jets) ){
    setCandidate(jets, closestToWMassIndices[0], lightQ_);
    match[TtSemiLepEvtPartons::LightQ] = closestToWMassIndices[0];
  }

  if( isValid(closestToWMassIndices[1], jets) ){
    setCandidate(jets, closestToWMassIndices[1], lightQBar_);
    match[TtSemiLepEvtPartons::LightQBar] = closestToWMassIndices[1];
  }

  for(unsigned idx=0; idx<maxPtIndices.size(); ++idx){
    // if this idx is not yet contained in the list of W mass candidates...
    if( std::find( closestToWMassIndices.begin(), closestToWMassIndices.end(), maxPtIndices[idx]) == closestToWMassIndices.end() ){
      // ...and if it is valid
      if( isValid(maxPtIndices[idx], jets) ){
	setCandidate(jets, maxPtIndices[idx], hadronicB_);
	match[TtSemiLepEvtPartons::HadB] = maxPtIndices[idx];
	break; // there should be no other cadidates!
      }
    }
  }

  if( isValid(lepB, jets) ){
    setCandidate(jets, lepB, leptonicB_);
    match[TtSemiLepEvtPartons::LepB] = lepB;
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

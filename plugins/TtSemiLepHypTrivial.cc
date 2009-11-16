#include "AnalysisDataFormats/TopObjects/interface/TtSemiLepEvtPartons.h"
#include "TopQuarkAnalysis/TopJetCombination/plugins/TtSemiLepHypTrivial.h"

TtSemiLepHypTrivial::TtSemiLepHypTrivial(const edm::ParameterSet& cfg):
  TtSemiLepHypothesis( cfg ),  
  maxNJets_(cfg.getParameter<int> ("maxNJets")),
  random_  (cfg.getParameter<bool>("random"  ))
{
  if(maxNJets_<4 && maxNJets_!=-1)
    throw cms::Exception("WrongConfig") 
      << "Parameter maxNJets can not be set to " << maxNJets_ << ". \n"
      << "It has to be larger than 4 or can be set to -1 to take all jets.";

  randNumGen_ = new TRandom3(0); // seed=0 -> computed via a TUUID object -> "unique in space and time"
}

TtSemiLepHypTrivial::~TtSemiLepHypTrivial()
{

  delete randNumGen_;

}

void
TtSemiLepHypTrivial::buildHypo(edm::Event& evt,
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
  for(unsigned int i=0; i<5; ++i)
    match.push_back(-1);

  std::vector<int> jetIndices;
  for(unsigned int i=0; i<maxNJets; ++i)
    jetIndices.push_back( i );
 
  int lightQ    = jetIndices[ randNumGen_->Integer( jetIndices.size() ) ];
  jetIndices.erase( std::find(jetIndices.begin(), jetIndices.end(), lightQ) );
 
  int lightQBar = jetIndices[ randNumGen_->Integer( jetIndices.size() ) ];
  jetIndices.erase( std::find(jetIndices.begin(), jetIndices.end(), lightQBar) );
 
  int hadB      = jetIndices[ randNumGen_->Integer( jetIndices.size() ) ];
  jetIndices.erase( std::find(jetIndices.begin(), jetIndices.end(), hadB) );
 
  int lepB      = jetIndices[0];
  if( jetIndices.size() > 1)
    lepB        = jetIndices[ randNumGen_->Integer( jetIndices.size() ) ];

  if(lightQ > lightQBar) {
    int itmp = lightQ;
    lightQ = lightQBar;
    lightQBar = itmp;
  }
  
  // -----------------------------------------------------
  // add jets
  // -----------------------------------------------------
  if( isValid(lightQ, jets) ){
    setCandidate(jets, lightQ, lightQ_, jetCorrectionLevel("wQuarkMix"));
    match[TtSemiLepEvtPartons::LightQ] = lightQ;
  }

  if( isValid(lightQBar, jets) ){
    setCandidate(jets, lightQBar, lightQBar_, jetCorrectionLevel("wQuarkMix"));
    match[TtSemiLepEvtPartons::LightQBar] = lightQBar;
  }

  if( isValid(hadB, jets) ){
    setCandidate(jets, hadB, hadronicB_, jetCorrectionLevel("bQuark"));
    match[TtSemiLepEvtPartons::HadB] = hadB;
  }
  
  if( isValid(lepB, jets) ){
    setCandidate(jets, lepB, leptonicB_, jetCorrectionLevel("bQuark"));
    match[TtSemiLepEvtPartons::LepB] = lepB;
  }

  // -----------------------------------------------------
  // add lepton
  // -----------------------------------------------------
  setCandidate(leps, 0, lepton_);
  match[TtSemiLepEvtPartons::Lepton] = 0;

  // -----------------------------------------------------
  // add neutrino
  // -----------------------------------------------------
  setCandidate(mets, 0, neutrino_);
}

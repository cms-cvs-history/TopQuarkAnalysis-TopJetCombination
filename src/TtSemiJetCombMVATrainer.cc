#include "TMath.h"

#include "PhysicsTools/MVATrainer/interface/HelperMacros.h"
#include "PhysicsTools/JetMCUtils/interface/combination.h"

#include "TopQuarkAnalysis/TopJetCombination/interface/TtSemiJetCombMVATrainer.h"
#include "TopQuarkAnalysis/TopJetCombination/interface/TtSemiJetCombEval.h"

TtSemiJetCombMVATrainer::TtSemiJetCombMVATrainer(const edm::ParameterSet& cfg):
  muons_   (cfg.getParameter<edm::InputTag>("muons")),
  jets_    (cfg.getParameter<edm::InputTag>("jets")),
  matching_(cfg.getParameter<edm::InputTag>("matching")),
  nJetsMax_(cfg.getParameter<int>("nJetsMax"))
{
}

TtSemiJetCombMVATrainer::~TtSemiJetCombMVATrainer()
{
}

void
TtSemiJetCombMVATrainer::analyze(const edm::Event& evt, const edm::EventSetup& setup)
{
  mvaComputer.update<TtSemiJetCombMVARcd>("trainer", setup, "ttSemiJetCombMVA");

  if(!mvaComputer) return;   // can occur in the last iteration when the MVATrainer is about to save the result

  edm::Handle<TopMuonCollection> muons; 
  evt.getByLabel(muons_, muons);

  edm::Handle<TopJetCollection> topJets;
  evt.getByLabel(jets_, topJets);

  edm::Handle< std::vector<int> > matching;
  evt.getByLabel(matching_, matching);

  for(unsigned int i = 0; i < matching->size(); ++i) {
    if( (*(matching))[i] < 0) return;  // skip events that were affected by the outlier rejection in the jet-parton matching
  }

  math::XYZTLorentzVector muon = (*(muons->begin())).p4();

  // analyze true and false jet combinations
  std::vector<int> jetIndices;
  for(unsigned int i=0; i<topJets->size(); ++i) {
    if(nJetsMax_>=4 && i==(unsigned int) nJetsMax_) break;
    jetIndices.push_back(i);
  }

  std::vector<int> combi;
  for(unsigned int i=0; i<matching->size(); ++i) 
    combi.push_back(i);

  do{
    for(int cnt=0; cnt<TMath::Factorial( matching->size() ); ++cnt){
      if(combi[0] < combi[1]) {  // take into account indistinguishability 
	                         // of the two jets from the hadr. W decay,
	                         // reduces combinatorics by a factor of 2
	TtSemiJetComb jetComb(*topJets, combi, muon);

	bool trueCombi = true;
	for(unsigned int i=0; i<matching->size(); ++i){
	  if(combi[i] != (*(matching))[i]) {
	    trueCombi = false;
	    break;
	  }
	}
	evaluateTtSemiJetComb(mvaComputer, jetComb, true, trueCombi);
      }
      next_permutation( combi.begin() , combi.end() );
    }    
  }
  while(stdcomb::next_combination( jetIndices.begin(), jetIndices.end(), combi.begin(), combi.end() ));
}

// implement the plugins for the trainer
// -> defines TtSemiJetCombMVAContainerSaveCondDB
// -> defines TtSemiJetCombMVASaveFile
// -> defines TtSemiJetCombMVATrainerLooper
MVA_TRAINER_IMPLEMENT(TtSemiJetCombMVA);
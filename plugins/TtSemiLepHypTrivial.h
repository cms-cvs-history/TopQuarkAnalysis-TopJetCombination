#ifndef TtSemiLepHypTrivial_h
#define TtSemiLepHypTrivial_h

#include "TopQuarkAnalysis/TopJetCombination/interface/TtSemiLepHypothesis.h"

#include "TRandom3.h"

class TtSemiLepHypTrivial : public TtSemiLepHypothesis  {

 public:

  explicit TtSemiLepHypTrivial(const edm::ParameterSet&);
  ~TtSemiLepHypTrivial();

 private:

  /// build the event hypothesis key
  virtual void buildKey() { key_= TtSemiLeptonicEvent::kTrivial; };  
  /// build event hypothesis from the reco objects of a semi-leptonic event 
  virtual void buildHypo(edm::Event&,
			 const edm::Handle<edm::View<reco::RecoCandidate> >&,
			 const edm::Handle<std::vector<pat::MET> >&,
			 const edm::Handle<std::vector<pat::Jet> >&,
			 std::vector<int>&, const unsigned int iComb);

 private:

  /// random number generator object
  TRandom3 *randNumGen_;

  int maxNJets_;
  bool random_;
};

#endif

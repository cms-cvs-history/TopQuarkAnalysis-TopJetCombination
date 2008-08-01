#ifndef TtSemiHypothesisKinFit_h
#define TtSemiHypothesisKinFit_h

#include "TopQuarkAnalysis/TopJetCombination/interface/TtSemiHypothesis.h"

class TtSemiHypothesisKinFit : public TtSemiHypothesis  {

 public:

  explicit TtSemiHypothesisKinFit(const edm::ParameterSet&);
  ~TtSemiHypothesisKinFit();

 private:
  
  /// build the event hypothesis key
  virtual void buildKey() { key_= TtSemiEvent::kKinFit; };  
  /// build event hypothesis from the reco objects of a semi-leptonic event 
  virtual void buildHypo(edm::Event&,
			 const edm::Handle<edm::View<reco::RecoCandidate> >&,
			 const edm::Handle<std::vector<pat::MET> >&,
			 const edm::Handle<std::vector<pat::Jet> >&,
			 std::vector<int>&);

  edm::InputTag status_;
  edm::InputTag hadB_;
  edm::InputTag hadP_;
  edm::InputTag hadQ_;
  edm::InputTag lepB_;
  edm::InputTag lepL_;
  edm::InputTag lepN_;

};

#endif

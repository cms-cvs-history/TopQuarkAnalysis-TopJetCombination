#ifndef TtSemiHypothesis_h
#define TtSemiHypothesis_h

#include <memory>
#include <vector>

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"

#include "DataFormats/Candidate/interface/CandidateWithRef.h"
#include "AnalysisDataFormats/TopObjects/interface/TtSemiEvent.h"
#include "DataFormats/Candidate/interface/ShallowCloneCandidate.h"
#include "DataFormats/Candidate/interface/NamedCompositeCandidate.h"


class TtSemiHypothesis : public edm::EDProducer {

 public:

  explicit TtSemiHypothesis(const edm::ParameterSet&);
  ~TtSemiHypothesis();

 protected:
  
  /// produce the event hypothesis as NamedCompositeCandidate and Key
  virtual void produce(edm::Event&, const edm::EventSetup&);
  /// return key
  int key() const { return key_; };
  /// return event hypothesis
  reco::NamedCompositeCandidate hypo();
  /// check if index is in valid range of selected jets
  bool isValid(const int& idx, const edm::Handle<std::vector<pat::Jet> >& jets){ return (0<=idx && idx<(int)jets->size()); };

  // -----------------------------------------
  // implemet the following two functions
  // for a concrete event hypothesis
  // -----------------------------------------

  /// build the event hypothesis key
  virtual void buildKey() = 0;
  /// build event hypothesis from the reco objects of a semi-leptonic event 
  virtual void buildHypo(const edm::Handle<edm::View<reco::RecoCandidate> >& lepton, 
			 const edm::Handle<std::vector<pat::MET> >& neutrino, 
			 const edm::Handle<std::vector<pat::Jet> >& jets, 
			 const edm::Handle<std::vector<int> >& jetPartonAssociation) = 0;

 protected:

  edm::InputTag jets_;
  edm::InputTag leps_;
  edm::InputTag mets_;
  edm::InputTag match_;  

  int key_;

  reco::ShallowCloneCandidate *lightQ_;
  reco::ShallowCloneCandidate *lightQBar_;
  reco::ShallowCloneCandidate *hadronicB_;
  reco::ShallowCloneCandidate *leptonicB_;
  reco::ShallowCloneCandidate *neutrino_;
  reco::ShallowCloneCandidate *lepton_;
};

#endif
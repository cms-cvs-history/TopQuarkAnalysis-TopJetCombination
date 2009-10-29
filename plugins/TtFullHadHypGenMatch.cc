#include "TopQuarkAnalysis/TopJetCombination/plugins/TtFullHadHypGenMatch.h"
#include "AnalysisDataFormats/TopObjects/interface/TtFullHadEvtPartons.h"
#include "DataFormats/Math/interface/deltaR.h"

TtFullHadHypGenMatch::TtFullHadHypGenMatch(const edm::ParameterSet& cfg):
  TtFullHadHypothesis( cfg ) 
{  
}

TtFullHadHypGenMatch::~TtFullHadHypGenMatch() { }

void
TtFullHadHypGenMatch::buildHypo(edm::Event& evt,
			        const edm::Handle<std::vector<pat::Jet> >& jets, 
			        std::vector<int>& match,
				const unsigned int iComb)
{
  // -----------------------------------------------------
  // add jets
  // -----------------------------------------------------
  for(unsigned idx=0; idx<match.size(); ++idx){
    if( isValid(match[idx], jets) ){
      switch(idx){
      case TtFullHadEvtPartons::LightQ:
	setCandidate(jets, match[idx], lightQ_   , jetCorrectionLevel_); break;
      case TtFullHadEvtPartons::LightQBar:
	setCandidate(jets, match[idx], lightQBar_, jetCorrectionLevel_); break;	
      case TtFullHadEvtPartons::B:
	setCandidate(jets, match[idx], b_   , jetCorrectionLevel_); break;
      case TtFullHadEvtPartons::LightP:
	setCandidate(jets, match[idx], lightP_, jetCorrectionLevel_); break;	
      case TtFullHadEvtPartons::LightPBar:
	setCandidate(jets, match[idx], lightPBar_   , jetCorrectionLevel_); break;
      case TtFullHadEvtPartons::BBar:
	setCandidate(jets, match[idx], bBar_, jetCorrectionLevel_); break;	
      }
    }
  }
}

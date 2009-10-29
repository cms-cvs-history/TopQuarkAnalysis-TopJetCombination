import FWCore.ParameterSet.Config as cms

#
# module to make the genMatch hypothesis
#
ttFullHadHypGenMatch = cms.EDProducer("TtFullHadHypGenMatch",
    jets  = cms.InputTag("selectedLayer1Jets"),    
    match = cms.InputTag("ttFullHadJetPartonMatch"), 
    jetCorrectionLevel = cms.string("abs")   
)



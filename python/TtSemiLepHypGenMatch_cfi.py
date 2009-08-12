import FWCore.ParameterSet.Config as cms

#
# module to make the genMatch hypothesis
#
ttSemiLepHypGenMatch = cms.EDProducer("TtSemiLepHypGenMatch",
    mets  = cms.InputTag("layer1METs"),
    leps  = cms.InputTag("selectedLayer1Muons"),
    jets  = cms.InputTag("selectedLayer1Jets"),
    match = cms.InputTag("ttSemiLepJetPartonMatch"),
    jetCorrectionLevel = cms.string("abs")
)



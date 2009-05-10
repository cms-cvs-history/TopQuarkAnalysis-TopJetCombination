import FWCore.ParameterSet.Config as cms

#
# module to make the mvaDiscriminator hypothesis
#
ttSemiLepHypMVADisc = cms.EDProducer("TtSemiLepHypMVADisc",
    mets  = cms.InputTag("layer1METs"),
    leps  = cms.InputTag("selectedLayer1Muons"),
    jets  = cms.InputTag("selectedLayer1Jets"),
    match = cms.InputTag("findTtSemiLepJetCombMVA")
)



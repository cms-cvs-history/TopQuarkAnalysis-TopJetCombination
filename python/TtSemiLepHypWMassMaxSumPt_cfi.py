import FWCore.ParameterSet.Config as cms

#
# module to make the wMassMaxSumPt hypothesis
#
ttSemiLepHypWMassMaxSumPt = cms.EDProducer("TtSemiLepHypWMassMaxSumPt",
    mets  = cms.InputTag("layer1METs"),
    leps  = cms.InputTag("selectedLayer1Muons"),
    jets  = cms.InputTag("selectedLayer1Jets"),
    maxNJets = cms.int32(4),
    wMass    = cms.double(80.413)
)



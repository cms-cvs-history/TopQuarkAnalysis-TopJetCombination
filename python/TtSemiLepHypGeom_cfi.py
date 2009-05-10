import FWCore.ParameterSet.Config as cms

#
# module to make the geom hypothesis
#
ttSemiLepHypGeom = cms.EDProducer("TtSemiLepHypGeom",
    mets  = cms.InputTag("layer1METs"),
    leps  = cms.InputTag("selectedLayer1Muons"),
    jets  = cms.InputTag("selectedLayer1Jets"),
    maxNJets  = cms.int32(4),
    useDeltaR = cms.bool(True)                               
)

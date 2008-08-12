import FWCore.ParameterSet.Config as cms

#
# module to make the geom hypothesis
#
ttSemiLepGeom = cms.EDProducer("TtSemiLepGeom",
    leps  = cms.InputTag("selectedLayer1Muons"),
    mets  = cms.InputTag("selectedLayer1METs"),
    jets  = cms.InputTag("selectedLayer1Jets"),
    maxNJets  = cms.int32(4),
    useDeltaR = cms.bool(True)                               
)

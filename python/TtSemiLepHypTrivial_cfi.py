import FWCore.ParameterSet.Config as cms

#
# module to make the trivial hypothesis
#
ttSemiLepHypTrivial = cms.EDProducer("TtSemiLepHypTrivial",
    ## met input 
    mets  = cms.InputTag("layer1METs"),
    ## jet input 
    jets  = cms.InputTag("selectedLayer1Jets"),
    ## lepton input
    leps  = cms.InputTag("selectedLayer1Muons"),
    ## maximal number of jets to be considered
    maxNJets  = cms.int32(4),
    ## randomly pick one of all possible jet combinations
    random = cms.bool(True)
)

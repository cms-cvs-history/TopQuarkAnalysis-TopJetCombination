import FWCore.ParameterSet.Config as cms

process = cms.Process("TQAF")
# initialize MessageLogger
process.load("FWCore.MessageLogger.MessageLogger_cfi")

# input file
process.load("TopQuarkAnalysis.Examples.test.RecoInput_cfi")

# TQAF ###
process.load("TopQuarkAnalysis.TopObjectProducers.tqafLayer1_full_cff")

process.load("TopQuarkAnalysis.TopEventProducers.sequences.ttGenEvent_cff")

process.load("TopQuarkAnalysis.TopEventProducers.producers.TtDecaySelection_cfi")

process.load("TopQuarkAnalysis.TopTools.TtSemiEvtJetPartonMatch_cfi")

process.load("TopQuarkAnalysis.TopJetCombination.TtSemiJetCombMVATrainer_Muons_cff")

process.p0 = cms.Path(process.tqafLayer1*process.makeGenEvt)
process.p1 = cms.Path(process.ttDecaySelection*process.ttSemiJetPartonMatch)
process.p2 = cms.Path(process.makeMVATraining)
process.MessageLogger.cerr.threshold = 'INFO'
process.ttDecaySelection.channel_1 = [0, 1, 0]


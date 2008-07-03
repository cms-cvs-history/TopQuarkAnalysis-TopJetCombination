import FWCore.ParameterSet.Config as cms

#-------------------------------------------------
# test cfg file for mva training for jet parton 
# association
#-------------------------------------------------
process = cms.Process("TEST")

## add message logger
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'

#-------------------------------------------------
# process configuration
#-------------------------------------------------

## define input
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:/afs/cern.ch/cms/PRS/top/cmssw-data/relval200-for-pat-testing/FullSimTTbar-210p5.1-AODSIM.100.root')
)

## define maximal number of events to loop over
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)

## configure process options
process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)

#-------------------------------------------------
# tqaf configuration
#-------------------------------------------------

## std sequence for tqaf layer1
process.load("TopQuarkAnalysis.TopObjectProducers.tqafLayer1_full_cff")

## std sequence for ttGenEvent
process.load("TopQuarkAnalysis.TopEventProducers.sequences.ttGenEvent_cff")

## configure ttDecaySelection
process.load("TopQuarkAnalysis.TopEventProducers.producers.TtDecaySelection_cfi")
process.ttDecaySelection.channel_1 = [0, 1, 0]

## configure jet parton matching
process.load("TopQuarkAnalysis.TopTools.TtSemiEvtJetPartonMatch_cfi")

## configure mva trainer
process.load("TopQuarkAnalysis.TopJetCombination.TtSemiJetCombMVATrainer_Muons_cff")

#-------------------------------------------------
# process paths;
#-------------------------------------------------

## make jet parton match
process.p0 = cms.Path(process.tqafLayer1 *
                      process.makeGenEvt *
                      process.ttDecaySelection *
                      process.ttSemiJetPartonMatch)

## make mva training
process.p1 = cms.Path(process.makeMVATraining)

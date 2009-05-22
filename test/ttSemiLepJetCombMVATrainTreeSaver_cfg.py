import FWCore.ParameterSet.Config as cms

process = cms.Process("TEST")

## add message logger
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.categories.append('TtSemiLepJetCombMVATrainer')
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
    default                    = cms.untracked.PSet( limit = cms.untracked.int32( 0) ),
    TtSemiLepJetCombMVATrainer = cms.untracked.PSet( limit = cms.untracked.int32(-1) )
)

## define input
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
    ## small test file on afs
    'file:/afs/cern.ch/cms/PRS/top/cmssw-data/relval200-for-pat-testing/FullSimTTBar-2_2_X_2008-11-03-STARTUP_V7-AODSIM.100.root'
    ## relval sample with 9,000 events
    #'/store/relval/CMSSW_2_2_6/RelValTTbar/GEN-SIM-RECO/IDEAL_V12_v1/0002/18B597E9-BB0B-DE11-8232-0030487A322E.root',
    #'/store/relval/CMSSW_2_2_6/RelValTTbar/GEN-SIM-RECO/IDEAL_V12_v1/0002/462BE681-490C-DE11-8558-000423D98750.root',
    #'/store/relval/CMSSW_2_2_6/RelValTTbar/GEN-SIM-RECO/IDEAL_V12_v1/0002/8646BC5A-C20B-DE11-96D8-001617DBD316.root',
    #'/store/relval/CMSSW_2_2_6/RelValTTbar/GEN-SIM-RECO/IDEAL_V12_v1/0002/BAA26726-C00B-DE11-858D-001D09F24EE3.root',
    #'/store/relval/CMSSW_2_2_6/RelValTTbar/GEN-SIM-RECO/IDEAL_V12_v1/0002/D4D3FE46-BD0B-DE11-8F6C-000423D9880C.root',
    #'/store/relval/CMSSW_2_2_6/RelValTTbar/GEN-SIM-RECO/IDEAL_V12_v1/0002/F802372C-BB0B-DE11-B812-000423D99AA2.root'
     ),
     skipEvents = cms.untracked.uint32(0)
)
## define maximal number of events to loop over
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)
## configure process options
process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(False)
)

## configure geometry & conditions
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string('STARTUP_V7::All')

## std sequence for pat
process.load("PhysicsTools.PatAlgos.patSequences_cff")

## std sequence for ttGenEvent
process.load("TopQuarkAnalysis.TopEventProducers.sequences.ttGenEvent_cff")

## configure ttGenEventFilters
process.load("TopQuarkAnalysis.TopEventProducers.sequences.ttGenEventFilters_cff")
process.ttSemiLeptonicFilter.allowedTopDecays.decayBranchA.electron = False
process.ttSemiLeptonicFilter.allowedTopDecays.decayBranchA.muon     = True
process.ttSemiLeptonicFilter.allowedTopDecays.decayBranchA.tau      = False

## configure jet-parton matching
process.load("TopQuarkAnalysis.TopTools.TtSemiLepJetPartonMatch_cfi")
#process.ttSemiLepJetPartonMatch.partonsToIgnore = ["LepB"]

## configure mva trainer
process.load("TopQuarkAnalysis.TopJetCombination.TtSemiLepJetCombMVATrainTreeSaver_cff")
## change maximum number of jets taken into account per event (default: 4)
#process.ttSemiLepJetPartonMatch .maxNJets = 5
#process.trainTtSemiLepJetCombMVA.maxNJets = process.ttSemiLepJetPartonMatch.maxNJets

## make trainer looper known to the process
from TopQuarkAnalysis.TopJetCombination.TtSemiLepJetCombMVATrainTreeSaver_cff import looper
process.looper = looper

process.p = cms.Path(process.makeGenEvt *
                     process.patDefaultSequence *
                     process.ttSemiLeptonicFilter *
                     process.ttSemiLepJetPartonMatch *
                     process.saveTtSemiLepJetCombMVATrainTree)
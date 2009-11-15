import FWCore.ParameterSet.Config as cms

#
# produce trivial hypothesis with all necessary 
# ingredients
#

## configure trivial hyothesis
from TopQuarkAnalysis.TopJetCombination.TtSemiLepHypTrivial_cfi import *

## make hypothesis
makeHypothesis_trivial = cms.Sequence(ttSemiLepHypTrivial)


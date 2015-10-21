import FWCore.ParameterSet.Config as cms

ExtraIsolationVariables = cms.PSet(
    fixedGridRhoFastjetAll                  = cms.InputTag("computeCorrectedIso", "fixedGridRhoFastjetAll"), # was RhoAll / kt6RhoAll
    fixedGridRhoFastjetAllCalo              = cms.InputTag("computeCorrectedIso", "fixedGridRhoFastjetAllCalo"), # was RhoAllCalo / kt6RhoAllCalo
    fixedGridRhoFastjetCentralChargedPileUp = cms.InputTag("computeCorrectedIso", "fixedGridRhoFastjetCentralChargedPileUp"), # was RhoPU / kt6RhoPU
    fixedGridRhoFastjetCentralNeutral       = cms.InputTag("computeCorrectedIso", "fixedGridRhoFastjetCentralNeutral"), # was RhoNeu05 / kt6RhoNeu05
    fixedGridRhoFastjetCentralCalo          = cms.InputTag("computeCorrectedIso", "fixedGridRhoFastjetCentralCalo"), # was RhoCentralCalo / kt6RhoCentralCalo
)

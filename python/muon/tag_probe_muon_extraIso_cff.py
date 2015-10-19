import FWCore.ParameterSet.Config as cms

ExtraIsolationVariables = cms.PSet(
    kt6RhoAll         = cms.InputTag("computeCorrectedIso", "fixedGridRhoFastjetAll"), # was RhoAll
    kt6RhoAllCalo     = cms.InputTag("computeCorrectedIso", "fixedGridRhoFastjetAllCalo"), # was RhoAllCalo
    kt6RhoPU          = cms.InputTag("computeCorrectedIso", "fixedGridRhoFastjetCentralChargedPileUp"), # was RhoPU
    kt6RhoNeu05       = cms.InputTag("computeCorrectedIso", "fixedGridRhoFastjetCentralNeutral"), # was RhoNeu05
    kt6RhoCentralCalo = cms.InputTag("computeCorrectedIso", "fixedGridRhoFastjetCentralCalo"), # was RhoCentralCalo
)

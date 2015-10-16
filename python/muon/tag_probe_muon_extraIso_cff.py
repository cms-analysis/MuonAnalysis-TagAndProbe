import FWCore.ParameterSet.Config as cms

ExtraIsolationVariables = cms.PSet(
    kt6RhoAll         = cms.InputTag("computeCorrectedIso", "RhoAll"), # fixedGridRhoFastjetAll
    kt6RhoAllCalo     = cms.InputTag("computeCorrectedIso", "RhoAllCalo"), # fixedGridRhoFastjetAllCalo
    kt6RhoPU          = cms.InputTag("computeCorrectedIso", "RhoPU"), # fixedGridRhoFastjetCentralChargedPileUp
    kt6RhoNeu05       = cms.InputTag("computeCorrectedIso", "RhoNeu05"), # fixedGridRhoFastjetCentralNeutral
    kt6RhoCentralCalo = cms.InputTag("computeCorrectedIso", "RhoCentralCalo"), # fixedGridRhoFastjetCentralCalo
    # kt6RhoNeu1        = cms.InputTag("computeCorrectedIso", "RhoNeu1"), # fixedGridRhoFastjetAll ### FIXME
)

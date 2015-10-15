import FWCore.ParameterSet.Config as cms

ExtraIsolationVariables = cms.PSet(
    kt6RhoAll         = cms.InputTag("computeCorrectedIso", "RhoAll"),
    kt6RhoAllCalo     = cms.InputTag("computeCorrectedIso", "RhoAllCalo"),
    kt6RhoPU          = cms.InputTag("computeCorrectedIso", "RhoPU"),
    kt6RhoNeu05       = cms.InputTag("computeCorrectedIso", "RhoNeu05"),
    kt6RhoNeu1        = cms.InputTag("computeCorrectedIso", "RhoNeu1"),
    fixedGridRhoAll = cms.InputTag("computeCorrectedIso", "fixedGridRhoAll"),
    fixedGridRhoFastjetAll = cms.InputTag("computeCorrectedIso", "fixedGridRhoFastjetAll"),
    fixedGridRhoFastjetAllCalo = cms.InputTag("computeCorrectedIso", "fixedGridRhoFastjetAllCalo"),
    fixedGridRhoFastjetCentralCalo = cms.InputTag("computeCorrectedIso", "fixedGridRhoFastjetCentralCalo"),
    fixedGridRhoFastjetCentralChargedPileUp = cms.InputTag("computeCorrectedIso", "fixedGridRhoFastjetCentralChargedPileUp"),
    fixedGridRhoFastjetCentralNeutral = cms.InputTag("computeCorrectedIso", "fixedGridRhoFastjetCentralNeutral")
)

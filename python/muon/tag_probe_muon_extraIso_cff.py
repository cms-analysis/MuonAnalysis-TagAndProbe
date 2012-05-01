import FWCore.ParameterSet.Config as cms

ExtraIsolationVariables = cms.PSet(
    kt6RhoAll         = cms.InputTag("computeCorrectedIso", "RhoAll"),
    kt6RhoAllCalo     = cms.InputTag("computeCorrectedIso", "RhoAllCalo"),
    kt6RhoPU          = cms.InputTag("computeCorrectedIso", "RhoPU"),
    kt6RhoNeu05       = cms.InputTag("computeCorrectedIso", "RhoNeu05"),
    kt6RhoNeu1        = cms.InputTag("computeCorrectedIso", "RhoNeu1"),
    ecalIsoRhoCorr    = cms.InputTag("computeCorrectedIso", "ecalIsoRhoCorrected"),
    hcalIsoRhoCorr    = cms.InputTag("computeCorrectedIso", "hcalIsoRhoCorrected"),
    combRelIsoRhoCorr = cms.InputTag("computeCorrectedIso", "combRelIsoRhoCorrected"),
    )

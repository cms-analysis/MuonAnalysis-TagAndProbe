import FWCore.ParameterSet.Config as cms

ExtraIsolationVariables = cms.PSet(
    kt6Rho            = cms.InputTag("computeCorrectedIso", "Rho"),
    ecalIsoRhoCorr    = cms.InputTag("computeCorrectedIso", "ecalIsoRhoCorrected"),
    hcalIsoRhoCorr    = cms.InputTag("computeCorrectedIso", "hcalIsoRhoCorrected"),
    combRelIsoRhoCorr = cms.InputTag("computeCorrectedIso", "combRelIsoRhoCorrected"),
    )

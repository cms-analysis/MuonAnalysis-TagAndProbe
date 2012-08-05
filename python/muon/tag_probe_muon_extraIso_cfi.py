import FWCore.ParameterSet.Config as cms

computeCorrectedIso = cms.EDProducer("ComputeIsoCorrections",
    probes = cms.InputTag("probeMuons"),
)

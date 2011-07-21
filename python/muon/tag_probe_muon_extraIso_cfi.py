import FWCore.ParameterSet.Config as cms

computeCorrectedIso = cms.EDProducer("ComputeIsoCorrections",
                                     probes = cms.InputTag("probeMuons"),
                                     EffAreaEcalBar = cms.double(0.074),
                                     EffAreaEcalEnd = cms.double(0.041),
                                     EffAreaHcalBar = cms.double(0.023),
                                     EffAreaHcalEnd = cms.double(0.032),
                                     )

import FWCore.ParameterSet.Config as cms

radialIso = cms.EDProducer("MuonRadialIso",
    probes = cms.InputTag("probeMuons"),
    pfCandidates = cms.InputTag("pfNoPileUpPFIso"),
    photonPtMin = cms.double(1.0),
    neutralHadPtMin = cms.double(1.0),
)

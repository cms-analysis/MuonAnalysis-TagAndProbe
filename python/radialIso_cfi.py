import FWCore.ParameterSet.Config as cms

# Santi's tune
radialIso = cms.EDProducer("MuonRadialIso",
    probes = cms.InputTag("probeMuons"),
    pfCandidates = cms.InputTag("pfNoPileUpMVAIso"),
    photonPtMin = cms.double(1.0),
    neutralHadPtMin = cms.double(1.0),
)

import FWCore.ParameterSet.Config as cms

# Si Xie's tune (but no overlap removal)
mvaIsoVariables = cms.EDProducer("MuonMVAIsoVariables",
    probes = cms.InputTag("probeMuons"),
    doOverlapRemoval = cms.bool(False), # remove electrons and muon taken from this lists from the sums
    goodElectrons = cms.InputTag("NO"), # 
    goodMuons     = cms.InputTag("NO"), # 
    pfCandidates = cms.InputTag("particleFlow"),
    vertices     = cms.InputTag("offlinePrimaryVertices"),
    dzCut = cms.double(0.2),
    photonPtMin = cms.double(-1),
    neutralHadPtMin = cms.double(-1),
)

import FWCore.ParameterSet.Config as cms

nverticesModule = cms.EDProducer("VertexMultiplicityCounter", 
    probes = cms.InputTag("tagMuons"),
    objects = cms.InputTag("offlinePrimaryVertices"),
    objectSelection = cms.string("!isFake && ndof > 4 && abs(z) <= 25 && position.Rho <= 2"),
)

njets15Module = cms.EDProducer("CandCleanedMultiplicityCounter", 
    pairs   = cms.InputTag("tpPairs"),
    objects = cms.InputTag("ak5PFJets"),
    objectSelection = cms.string("abs(eta) < 5 && pt > 15"), 
    minTagObjDR   = cms.double(0.3),
    minProbeObjDR = cms.double(0.3),
)

njets30Module = cms.EDProducer("CandCleanedMultiplicityCounter", 
    pairs   = cms.InputTag("tpPairs"),
    objects = cms.InputTag("ak5PFJets"),
    objectSelection = cms.string("abs(eta) < 5 && pt > 30"), 
    minTagObjDR   = cms.double(0.3),
    minProbeObjDR = cms.double(0.3),
)

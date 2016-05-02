import FWCore.ParameterSet.Config as cms

centralityInfo = cms.EDProducer("HeavyIonCentralityInfo",
	src = cms.InputTag("tagMuons"),
        CentralitySrc = cms.InputTag("hiCentrality"),
)

centralityBinInfo = cms.EDProducer("HeavyIonCentralityBinInfo",
	src = cms.InputTag("tagMuons"),
        CentralityBinSrc = cms.InputTag("centralityBin","HFtowers"),
)

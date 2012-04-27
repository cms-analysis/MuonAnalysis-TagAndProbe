import FWCore.ParameterSet.Config as cms

from MuonAnalysis.TagAndProbe.common_modules_cff import nverticesModule
nverticesModuleProbes = nverticesModule.clone(probes = 'probeMuons')

probeRecoMuons = cms.EDFilter("PATMuonRefSelector",
    src = cms.InputTag("probeMuons"),
    cut = cms.string("isGlobalMuon || isPFMuon || isTrackerMuon"),
)

diMuPairs = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("probeRecoMuons@+ probeRecoMuons@-"),
    cut   = cms.string("mass > 1"),
)


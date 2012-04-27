import FWCore.ParameterSet.Config as cms

from MuonAnalysis.TagAndProbe.fakerate_common_cff import *

zForFakes = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("tagMuons@+ probeRecoMuons@-"),
    cut   = cms.string("60 < mass < 140 && daughter(1).pt > 10 && "+
                       "daughter(1).masterClone.isGlobalMuon && daughter(1).masterClone.numberOfMatches > 1 &&"+
                       "daughter(0).masterClone.isolationR03.sumPt/daughter(0).pt + daughter(1).masterClone.isolationR03.sumPt/daughter(1).pt < 0.3"),
)

zPlusProbe = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("zForFakes probeRecoMuons"),
    cut   = cms.string("deltaR(daughter(0).daughter(0).eta, daughter(0).daughter(0).phi, daughter(1).eta, daughter(1).phi) > 0.3 && "+
                       "deltaR(daughter(0).daughter(1).eta, daughter(0).daughter(1).phi, daughter(1).eta, daughter(1).phi) > 0.3"),
    checkCharge = cms.bool(False),
)

zPlusProbeFilter = cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag("zPlusProbe"),
    minNumber = cms.uint32(1),
)

metVetoForZ = cms.EDFilter("CandViewSelector",
    src = cms.InputTag("pfMet"),
    cut = cms.string("pt > 25"),
    filter = cms.bool(True),
)
fourthLepVetoForZ = cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag("probeRecoMuons"),
    minNumber = cms.uint32(4),
)


zPlusProbeVetoSequence = cms.Sequence(~metVetoForZ + ~fourthLepVetoForZ)

zPlusProbeSequence = cms.Sequence(
    probeRecoMuons +
    zForFakes * zPlusProbe * zPlusProbeFilter +
    zPlusProbeVetoSequence
)


ZPlusProbeTagVariables = cms.PSet(
    z1Mass = cms.string("mass"),
    mu1Eta = cms.string("daughter(0).eta"),
    mu1Pt  = cms.string("daughter(0).pt"),
    mu2Eta = cms.string("daughter(1).eta"),
    mu2Pt  = cms.string("daughter(1).pt"),
)


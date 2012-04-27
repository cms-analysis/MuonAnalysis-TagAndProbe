import FWCore.ParameterSet.Config as cms

from MuonAnalysis.TagAndProbe.fakerate_common_cff import *

wMuNu = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("tagMuons pfMet"),
    cut   = cms.string("daughter(1).pt > 20 && sqrt(2*daughter(0).pt*daughter(1).pt*(1 - cos(daughter(0).phi - daughter(1).phi))) > 25"),
    checkCharge = cms.bool(False),
)

wPlusProbe = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("wMuNu probeRecoMuons"),
    cut   = cms.string("deltaR(daughter(0).daughter(0).eta, daughter(0).daughter(0).phi, daughter(1).eta, daughter(1).phi) > 0.5"),
    checkCharge = cms.bool(False),
)

wPlusProbeFilter = cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag("wPlusProbe"),
    minNumber = cms.uint32(1),
)

zVetoForW = cms.EDFilter("CandViewSelector",
    src = cms.InputTag("diMuPairs"),
    cut = cms.string("(mass < 12 || abs(mass-91.1876) < 30)"),
    filter = cms.bool(True),
)

jetsForVeto = cms.EDFilter("PFJetSelector",
    src = cms.InputTag("ak5PFJets"),
    cut = cms.string("pt > 30"),
)
jetVeto = cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag("jetsForVeto"),
    minNumber = cms.uint32(2),
)


wPlusProbeVetoSequence = cms.Sequence(~zVetoForW + jetsForVeto * ~jetVeto)

wPlusProbeSequence = cms.Sequence(
    probeRecoMuons * diMuPairs +
    wMuNu * 
    wPlusProbe * wPlusProbeFilter +
    wPlusProbeVetoSequence
)


WPlusProbeTagVariables = cms.PSet(
    pfMET = cms.string("daughter(1).pt"),
    muEta = cms.string("daughter(0).eta"),
    muPt  = cms.string("daughter(0).pt"),
    WMT   = cms.string("sqrt(2*daughter(0).pt*daughter(1).pt*(1 - cos(daughter(0).phi - daughter(1).phi)))"),
    trkRelIso = cms.string("daughter(0).masterClone.isolationR03.sumPt/daughter(0).pt"),
)


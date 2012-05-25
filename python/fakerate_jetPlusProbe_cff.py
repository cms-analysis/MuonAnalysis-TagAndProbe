import FWCore.ParameterSet.Config as cms

from MuonAnalysis.TagAndProbe.fakerate_common_cff import *

realJets = cms.EDFilter("CandViewSelector",
        src = cms.InputTag("ak5PFJets"),
        cut = cms.string("pt > 15 && (chargedMultiplicity+neutralMultiplicity) > 1 && chargedMuEnergyFraction < 0.25"),
        filter = cms.bool(True),
)
leadingJet = cms.EDFilter("LargestPtCandViewSelector",
    src = cms.InputTag("realJets"), 
    maxNumber = cms.uint32(1)
)

jetPlusProbe = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("leadingJet probeRecoMuons"),
    cut   = cms.string("deltaR(daughter(0).eta, daughter(0).phi, daughter(1).eta, daughter(1).phi) > 1.0 && daughter(1).pt > 10"),
    checkCharge = cms.bool(False)
)

jetPlusProbeFilter = cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag("jetPlusProbe"),
    minNumber = cms.uint32(1),
)

diLepVetoForJets = cms.EDFilter("CandViewSelector",
    src = cms.InputTag("diMuPairs"),
    cut = cms.string("(mass < 12 || abs(mass-91.1876) < 30)"),
    filter = cms.bool(True),
)
wMuNuForVeto = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("probeMuons pfMet"),
    cut   = cms.string("daughter(1).pt > 20 || sqrt(2*daughter(0).pt*daughter(1).pt*(1 - cos(daughter(0).phi - daughter(1).phi))) > 20"),
    checkCharge = cms.bool(False),
)
wMuNuVeto = cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag("wMuNuForVeto"),
    minNumber = cms.uint32(1),
)

jetPlusProbeVetoSequence = cms.Sequence(~diLepVetoForJets + wMuNuForVeto * ~wMuNuVeto)

jetPlusProbeSequence = cms.Sequence(
    probeRecoMuons * diMuPairs +
    jetPlusProbeVetoSequence   +
    realJets * leadingJet * jetPlusProbe * jetPlusProbeFilter
)


JetPlusProbeTagVariables = cms.PSet(
    pt  = cms.string("pt"),
    eta = cms.string("eta"),
)


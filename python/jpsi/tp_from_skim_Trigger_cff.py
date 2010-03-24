import FWCore.ParameterSet.Config as cms

from MuonAnalysis.TagAndProbe.jpsi.tp_from_skim_common_cff import *

anyProbeMuons = cms.EDFilter("PATMuonRefSelector",
    src = cms.InputTag("patMuons"),
    cut = cms.string("track.isNonnull && "+TRACK_CUTS+" && "+PT_ETA_CUTS), # pick everything here, select Glb/Trk later
)

tpGlbAny = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("tagMuons1Mu@+ anyProbeMuons@-"), # charge coniugate states are implied
    cut   = cms.string("%f < mass < %f" % MASS_RANGE),
)

histoTrigger = tnpTreeProducer.clone(
    tagProbePairs = cms.InputTag("tpGlbAny"),
    # choice of what defines a 'passing' probe
    flags = cms.PSet(
        Glb = cms.string(PASSING_GLB_CUT),
        Trk = cms.string(PASSING_TM_CUT),
        HLTMu3     = cms.string("!triggerObjectMatchesByFilter('hltSingleMu3L3Filtered3').empty()"),
        L1DiMuOpen = cms.string("!triggerObjectMatchesByFilter('hltDoubleMuLevel1PathL1OpenFiltered').empty()"),
    ),
    probeMatches  = cms.InputTag("muMcMatch"),
    allProbes = cms.InputTag("anyProbeMuons"),
)

tnpSequenceTrigger = cms.Sequence(
    anyProbeMuons *
    tpGlbAny *
    histoTrigger
)

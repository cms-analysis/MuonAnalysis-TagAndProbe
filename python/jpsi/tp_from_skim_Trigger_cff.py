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
        ## Selections
        Cal = cms.string("isCaloMuon"),
        Glb = cms.string(PASSING_GLB_CUT),
        TM  = cms.string(PASSING_TMI_CUT),
        POG_Glb    = cms.string("isGlobalMuon"),
        POG_GlbPT  = cms.string("isGlobalMuon  && muonID('GlobalMuonPromptTight')"),
        POG_TMA    = cms.string("isTrackerMuon && muonID('TrackerMuonArbitrated')"),
        POG_TMLSAT = cms.string("isTrackerMuon && muonID('TMLastStationAngTight')"),
        VBTFLike   = cms.string(PASSING_VBTFLIKE_CUTS),
        ## Triggers
        L1MuOpen   = cms.string("!triggerObjectMatchesByFilter('hltL1MuOpenL1Filtered0').empty()"),
        L2Mu0      = cms.string("!triggerObjectMatchesByFilter('hltL2Mu0L2Filtered0').empty()"),
        L2Mu3      = cms.string("!triggerObjectMatchesByFilter('hltSingleMu3L2Filtered3').empty()"),
        Mu3        = cms.string("!triggerObjectMatchesByFilter('hltSingleMu3L3Filtered3').empty()"),
        Mu5        = cms.string("!triggerObjectMatchesByFilter('hltSingleMu5L3Filtered5').empty()"),
        DoubleMu0  = cms.string("!triggerObjectMatchesByFilter('hltDiMuonL3PreFiltered').empty()"),
        DoubleMu3  = cms.string("!triggerObjectMatchesByFilter('hltDiMuonL3PreFiltered0').empty()"),
        L1DoubleMuOpen = cms.string("!triggerObjectMatchesByFilter('hltDoubleMuLevel1PathL1OpenFiltered').empty()"),
        ## Acceptance
        Acc_JPsi   = cms.string(JPSI_ACCEPTANCE_CUT),
    ),
    probeMatches  = cms.InputTag("muMcMatch"),
    allProbes = cms.InputTag("anyProbeMuons"),
)
histoTrigger.variables.caloCompatibility = cms.string("caloCompatibility")

tnpSequenceTrigger = cms.Sequence(
    anyProbeMuons *
    tpGlbAny *
    histoTrigger
)

def Add_CSCTF_Flags(histoTrigger):
    histoTrigger.flags.Ext_L1_Any   =  cms.string("userInt('muonL1MatchExtended') > 0")
    histoTrigger.flags.Ext_L1_Geom  =  cms.string("userInt('muonL1MatchExtended') >= 10")
    histoTrigger.flags.Ext_L1_CSC   =  cms.string("userInt('muonL1MatchExtended') != 0 && userInt('muonL1MatchExtended') != 10")
    histoTrigger.flags.Ext_L1_CSC_S =  cms.string("userInt('muonL1MatchExtended:cscMode') == 11")
    histoTrigger.flags.Ext_L1_CSC_M =  cms.string("userInt('muonL1MatchExtended') == 4 || userInt('muonL1MatchExtended') == 14")

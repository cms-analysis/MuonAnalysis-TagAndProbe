import FWCore.ParameterSet.Config as cms

from MuonAnalysis.TagAndProbe.jpsi.tp_from_skim_common_cff import *

##    ____                   _____               _      ____            _               
##   | __ )  __ _ _ __ ___  |_   _| __ __ _  ___| | __ |  _ \ _ __ ___ | |__   ___  ___ 
##   |  _ \ / _` | '__/ _ \   | || '__/ _` |/ __| |/ / | |_) | '__/ _ \| '_ \ / _ \/ __|
##   | |_) | (_| | | |  __/   | || | | (_| | (__|   <  |  __/| | | (_) | |_) |  __/\__ \
##   |____/ \__,_|_|  \___|   |_||_|  \__,_|\___|_|\_\ |_|   |_|  \___/|_.__/ \___||___/
##                                                                                      
##   
tkProbes = cms.EDFilter("CandViewRefSelector",
    src = cms.InputTag("tkTracks"),
    cut = cms.string(PT_ETA_CUTS),
)

##     ____      _       __  __                     ____            _               
##    / ___|__ _| | ___ |  \/  |_   _  ___  _ __   |  _ \ _ __ ___ | |__   ___  ___ 
##   | |   / _` | |/ _ \| |\/| | | | |/ _ \| '_ \  | |_) | '__/ _ \| '_ \ / _ \/ __|
##   | |__| (_| | | (_) | |  | | |_| | (_) | | | | |  __/| | | (_) | |_) |  __/\__ \
##    \____\__,_|_|\___/|_|  |_|\__,_|\___/|_| |_| |_|   |_|  \___/|_.__/ \___||___/
##                                                                                  
##   
calProbes = cms.EDFilter("PATMuonRefSelector",
    src = cms.InputTag("patMuons"),
    cut = cms.string("track.isNonnull && isCaloMuon && "+PT_ETA_CUTS+" && "+TRACK_CUTS),
)

allProbesMuonID = cms.Sequence(
    tkProbes +
    calProbes 
)

##    ____               _               _____               _      ____            _               
##   |  _ \ __ _ ___ ___(_)_ __   __ _  |_   _| __ __ _  ___| | __ |  _ \ _ __ ___ | |__   ___  ___ 
##   | |_) / _` / __/ __| | '_ \ / _` |   | || '__/ _` |/ __| |/ / | |_) | '__/ _ \| '_ \ / _ \/ __|
##   |  __/ (_| \__ \__ \ | | | | (_| |   | || | | (_| | (__|   <  |  __/| | | (_) | |_) |  __/\__ \
##   |_|   \__,_|___/___/_|_| |_|\__, |   |_||_|  \__,_|\___|_|\_\ |_|   |_|  \___/|_.__/ \___||___/
##                               |___/                                                              
##   

## As we did embed the tracks in the pat::Muons, and skim the generalTracks
## We can't do equality by reference, and we need a matching
muonsGlb = cms.EDFilter("PATMuonRefSelector",
    src = cms.InputTag("patMuons"),
    cut = cms.string(PASSING_GLB_CUT),
)
muonsTM = muonsGlb.clone(cut = PASSING_TM_CUT)

muonsPOGGlb    = muonsGlb.clone(cut = "isGlobalMuon")
muonsPOGGlbPT  = muonsGlb.clone(cut = "isGlobalMuon && muonID('GlobalMuonPromptTight')")
muonsPOGTMA    = muonsGlb.clone(cut = "isTrackerMuon && muonID('TrackerMuonArbitrated')")
muonsPOGTMLSAT = muonsGlb.clone(cut = "isTrackerMuon && muonID('TMLastStationAngTight')")
muonsVBTFLike  = muonsGlb.clone(cut = PASSING_VBTFLIKE_CUTS)

tkToGlbMatch = cms.EDProducer("MatcherUsingTracks",
    src     = cms.InputTag("tkTracks"), # all tracks are available for matching
    matched = cms.InputTag("muonsGlb"), # to all global muons
    algorithm = cms.string("byDirectComparison"), # check that they
    srcTrack = cms.string("tracker"),             # have the same 
    srcState = cms.string("atVertex"),            # tracker track
    matchedTrack = cms.string("tracker"),         # can't check ref
    matchedState = cms.string("atVertex"),        # because of the
    maxDeltaR        = cms.double(0.01),          # embedding.
    maxDeltaLocalPos = cms.double(0.01),
    maxDeltaPtRel    = cms.double(0.01),
    sortBy           = cms.string("deltaR"),
)
tkPassingGlb = cms.EDProducer("MatchedCandidateSelector",
    src   = cms.InputTag("tkProbes"),
    match = cms.InputTag("tkToGlbMatch"),
)
tkToTMMatch = tkToGlbMatch.clone(matched = 'muonsTM')
tkPassingTM = tkPassingGlb.clone( match = 'tkToTMMatch')

tkToStaMatch = tkToGlbMatch.clone(
    maxDeltaR        = cms.double(1.),   # large range in DR
    maxDeltaEta      = cms.double(0.4),  # small in eta, which is more precise
    maxDeltaLocalPos = cms.double(100),
    maxDeltaPtRel    = cms.double(5),
)

tkToPOGGlbMatch    = tkToGlbMatch.clone(matched = 'muonsPOGGlb')
tkToPOGGlbPTMatch  = tkToGlbMatch.clone(matched = 'muonsPOGGlbPT')
tkToPOGTMAMatch    = tkToGlbMatch.clone(matched = 'muonsPOGTMA')
tkToPOGTMLSATMatch = tkToGlbMatch.clone(matched = 'muonsPOGTMLSAT')
tkToVBTFLikeMatch  = tkToGlbMatch.clone(matched = 'muonsVBTFLike')
tkToPOGStaMatch    = tkToStaMatch.clone(matched = 'staTracks')
tkToPOGStaVHMatch  = tkToStaMatch.clone(matched = 'staTracksValidHits')
tkPassingPOGGlb    = tkPassingGlb.clone(match = 'tkToPOGGlbMatch')
tkPassingPOGGlbPT  = tkPassingGlb.clone(match = 'tkToPOGGlbPTMatch')
tkPassingPOGSta    = tkPassingGlb.clone(match = 'tkToPOGStaMatch')
tkPassingPOGStaVH  = tkPassingGlb.clone(match = 'tkToPOGStaVHMatch')
tkPassingPOGTMA    = tkPassingGlb.clone(match = 'tkToPOGTMAMatch')
tkPassingPOGTMLSAT = tkPassingGlb.clone(match = 'tkToPOGTMLSATMatch')
tkPassingVBTFLike  = tkPassingGlb.clone(match = 'tkToVBTFLikeMatch')

tkPassingProbes = cms.Sequence(
    muonsGlb * tkToGlbMatch * tkPassingGlb +
    muonsTM  * tkToTMMatch  * tkPassingTM 
)
tkPassingProbesPOG = cms.Sequence(
    muonsPOGGlb    * tkToPOGGlbMatch    * tkPassingPOGGlb    +
    muonsPOGGlbPT  * tkToPOGGlbPTMatch  * tkPassingPOGGlbPT  +
    muonsPOGTMA    * tkToPOGTMAMatch    * tkPassingPOGTMA    +
    muonsPOGTMLSAT * tkToPOGTMLSATMatch * tkPassingPOGTMLSAT +
                     tkToPOGStaMatch    * tkPassingPOGSta    +
                     tkToPOGStaVHMatch  * tkPassingPOGStaVH  +
    muonsVBTFLike  * tkToVBTFLikeMatch  * tkPassingVBTFLike
)


allPassingProbesMuonID = cms.Sequence(
    tkPassingProbes    +
    tkPassingProbesPOG
)
##    _____ ___   ____    ____       _          
##   |_   _( _ ) |  _ \  |  _ \ __ _(_)_ __ ___ 
##     | | / _ \/\ |_) | | |_) / _` | | '__/ __|
##     | || (_>  <  __/  |  __/ (_| | | |  \__ \
##     |_| \___/\/_|     |_|   \__,_|_|_|  |___/
##                                              
##   
tpGlbTk = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("tagMuons1Mu@+ tkProbes@-"), # charge coniugate states are implied
    cut   = cms.string("%f < mass < %f" % MASS_RANGE),
)
tpGlbTkVtxFit = cms.EDProducer("KalmanVertexFitCompositeCandProducer",
    src = cms.InputTag("tpGlbTk"),
)
tpGlbTkVtxCut = cms.EDFilter("CandViewRefSelector",
    src = cms.InputTag("tpGlbTkVtxFit"),
    cut = cms.string("chi2prob(vertexChi2, vertexNdof) > 0.01"), # 0.002 = chi2 cut at 10
)


tpGlbCal = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("tagMuons1Mu@+ calProbes@-"), # charge coniugate states are implied
    cut   = cms.string("%f < mass < %f" % MASS_RANGE),
)
tpGlbCalVtxFit = tpGlbTkVtxFit.clone(src = 'tpGlbCal')
tpGlbCalVtxCut = tpGlbTkVtxFit.clone(src = 'tpGlbCalVtxFit')

allTPPairsMuonID = cms.Sequence(
    tpGlbTk + tpGlbCal
)

##    __  __  ____   __  __       _       _               
##   |  \/  |/ ___| |  \/  | __ _| |_ ___| |__   ___  ___ 
##   | |\/| | |     | |\/| |/ _` | __/ __| '_ \ / _ \/ __|
##   | |  | | |___  | |  | | (_| | || (__| | | |  __/\__ \
##   |_|  |_|\____| |_|  |_|\__,_|\__\___|_| |_|\___||___/
##   
tkMcMatch  = muMcMatch.clone(src = "tkTracks")
allMcMatchesMuonID = cms.Sequence(tkMcMatch)


##    _____           _       _ ____            _            _   _ _____            _      
##   |_   _|_ _  __ _( )_ __ ( )  _ \ _ __ ___ | |__   ___  | \ | |_   _|   _ _ __ | | ___ 
##     | |/ _` |/ _` |/| '_ \|/| |_) | '__/ _ \| '_ \ / _ \ |  \| | | || | | | '_ \| |/ _ \
##     | | (_| | (_| | | | | | |  __/| | | (_) | |_) |  __/ | |\  | | || |_| | |_) | |  __/
##     |_|\__,_|\__, | |_| |_| |_|   |_|  \___/|_.__/ \___| |_| \_| |_| \__,_| .__/|_|\___|
##              |___/                                                        |_|           
##
#####
## Mu from Tk
histoMuFromTk = tnpTreeProducer.clone(
    tagProbePairs = cms.InputTag("tpGlbTk"),
    # choice of what defines a 'passing' probe
    flags = cms.PSet(
        Glb = cms.InputTag("tkPassingGlb"),
        TM  = cms.InputTag("tkPassingTM"),
        POG_Glb    = cms.InputTag("tkPassingPOGGlb"),
        POG_GlbPT  = cms.InputTag("tkPassingPOGGlbPT"),
        POG_Sta    = cms.InputTag("tkPassingPOGSta"),
        POG_StaVH  = cms.InputTag("tkPassingPOGStaVH"),
        POG_TMA    = cms.InputTag("tkPassingPOGTMA"),
        POG_TMLSAT = cms.InputTag("tkPassingPOGTMLSAT"),
        VBTFLike   = cms.InputTag("tkPassingVBTFLike"),
        Acc_JPsi   = cms.string(JPSI_ACCEPTANCE_CUT),
    ),
    ## These two MC things depend on the specific choice of probes
    probeMatches  = cms.InputTag("tkMcMatch"),
    allProbes     = cms.InputTag("tkProbes"), # NO 'unbias' efficiency on skims
)

#####
## Mu from Cal
histoMuFromCal = tnpTreeProducer.clone(
    tagProbePairs = cms.InputTag("tpGlbCal"),
    # choice of what defines a 'passing' probe
    flags = cms.PSet(
        Glb = cms.string(PASSING_GLB_CUT),
        TM  = cms.string(PASSING_TM_CUT),
        POG_Glb    = cms.string("isGlobalMuon"),
        POG_GlbPT  = cms.string("isGlobalMuon  && muonID('GlobalMuonPromptTight')"),
        POG_TMA    = cms.string("isTrackerMuon && muonID('TrackerMuonArbitrated')"),
        POG_TMLSAT = cms.string("isTrackerMuon && muonID('TMLastStationAngTight')"),
        VBTFLike   = cms.string(PASSING_VBTFLIKE_CUTS),
        Acc_JPsi   = cms.string(JPSI_ACCEPTANCE_CUT),
    ),
    ## These two MC things depend on the specific choice of probes
    probeMatches  = cms.InputTag("muMcMatch"),
    allProbes     = cms.InputTag("calProbes"),   
)

#####
## Mu from Tk
histoMuFromTkVtx = tnpTreeProducer.clone(
    tagProbePairs = cms.InputTag("tpGlbTkVtxCut"),
    # choice of what defines a 'passing' probe
    flags = cms.PSet(
        Glb = cms.InputTag("tkPassingGlb"),
        TM  = cms.InputTag("tkPassingTM"),
        POG_Glb    = cms.InputTag("tkPassingPOGGlb"),
        POG_GlbPT  = cms.InputTag("tkPassingPOGGlbPT"),
        POG_Sta    = cms.InputTag("tkPassingPOGSta"),
        POG_StaVH  = cms.InputTag("tkPassingPOGStaVH"),
        POG_TMA    = cms.InputTag("tkPassingPOGTMA"),
        POG_TMLSAT = cms.InputTag("tkPassingPOGTMLSAT"),
        VBTFLike   = cms.InputTag("tkPassingVBTFLike"),
        Acc_JPsi   = cms.string(JPSI_ACCEPTANCE_CUT),
    ),
    # chi2 of constrained fit (-1 if it failed)
    pairVariables = cms.PSet(vtxChi2 = cms.string("vertexChi2")),
    pairFlags     = cms.PSet(),
    ## These two MC things depend on the specific choice of probes
    probeMatches  = cms.InputTag("tkMcMatch"),
    allProbes     = cms.InputTag("tkProbes"), # NO 'unbias' efficiency on skims
)
histoMuFromCal.variables.caloCompatibility = cms.string("caloCompatibility")

#####
## Mu from Cal
histoMuFromCalVtx = tnpTreeProducer.clone(
    tagProbePairs = cms.InputTag("tpGlbCalVtxCut"),
    # choice of what defines a 'passing' probe
    flags = cms.PSet(
        Glb = cms.string(PASSING_GLB_CUT),
        TM  = cms.string(PASSING_TM_CUT),
        POG_Glb    = cms.string("isGlobalMuon"),
        POG_GlbPT  = cms.string("isGlobalMuon  && muonID('GlobalMuonPromptTight')"),
        POG_TMA    = cms.string("isTrackerMuon && muonID('TrackerMuonArbitrated')"),
        POG_TMLSAT = cms.string("isTrackerMuon && muonID('TMLastStationAngTight')"),
        VBTFLike   = cms.string(PASSING_VBTFLIKE_CUTS),
        Acc_JPsi   = cms.string(JPSI_ACCEPTANCE_CUT),
    ),
    # chi2 of constrained fit (-1 if it failed)
    pairVariables = cms.PSet(vtxChi2 = cms.string("vertexChi2")),
    pairFlags     = cms.PSet(),
    ## These two MC things depend on the specific choice of probes
    probeMatches  = cms.InputTag("muMcMatch"),
    allProbes     = cms.InputTag("calProbes"),
)
histoMuFromCalVtx.variables.caloCompatibility = cms.string("caloCompatibility")

allTPHistosMuonID = cms.Sequence(
    histoMuFromTk  +
    histoMuFromCal 
)

tnpSequenceMuonID = cms.Sequence(
    allProbesMuonID *
    allPassingProbesMuonID *
    allTPPairsMuonID   *
    allMcMatchesMuonID * 
    allTPHistosMuonID
)

tnpSequenceMuonIDVtx = cms.Sequence(
    tpGlbTkVtxFit  + tpGlbTkVtxCut  +
    tpGlbCalVtxFit + tpGlbCalVtxCut +
    histoMuFromTkVtx  +
    histoMuFromCalVtx 
)




import FWCore.ParameterSet.Config as cms

from MuonAnalysis.TagAndProbe.jpsi.tp_from_skim_common_cff import *

##    ____  _                  _    _    _                    ____            _               
##   / ___|| |_ __ _ _ __   __| |  / \  | | ___  _ __   ___  |  _ \ _ __ ___ | |__   ___  ___ 
##   \___ \| __/ _` | '_ \ / _` | / _ \ | |/ _ \| '_ \ / _ \ | |_) | '__/ _ \| '_ \ / _ \/ __|
##    ___) | || (_| | | | | (_| |/ ___ \| | (_) | | | |  __/ |  __/| | | (_) | |_) |  __/\__ \
##   |____/ \__\__,_|_| |_|\__,_/_/   \_\_|\___/|_| |_|\___| |_|   |_|  \___/|_.__/ \___||___/
##                                                                                            
##   
staProbes = cms.EDProducer("CandViewRefSelector",
    src = cms.InputTag("staTracks"),
    cut = cms.string(PT_ETA_CUTS),
)

tkTracksNoJPsi = cms.EDProducer("CandidateResonanceInefficiencyCreator",
    src = cms.InputTag("tkTracks"),
    tags = cms.InputTag("tagMuons2Mu"),
    mass    = cms.double(3.096),
    massMin = cms.double(2.85), ## Should cut away
    massMax = cms.double(3.25), ## 99.5% of signal
    onlyBestMatch = cms.bool(False),
    outputMode = cms.string("RefToBaseVector"),
)
tkTracksNoBestJPsi = tkTracksNoJPsi.clone(onlyBestMatch = True)

##    ____               _               ____            _                   _____               _    _             
##   |  _ \ __ _ ___ ___(_)_ __   __ _  |  _ \ _ __ ___ | |__   ___  ___ _  |_   _| __ __ _  ___| | _(_)_ __   __ _ 
##   | |_) / _` / __/ __| | '_ \ / _` | | |_) | '__/ _ \| '_ \ / _ \/ __(_)   | || '__/ _` |/ __| |/ / | '_ \ / _` |
##   |  __/ (_| \__ \__ \ | | | | (_| | |  __/| | | (_) | |_) |  __/\__ \_    | || | | (_| | (__|   <| | | | | (_| |
##   |_|   \__,_|___/___/_|_| |_|\__, | |_|   |_|  \___/|_.__/ \___||___(_)   |_||_|  \__,_|\___|_|\_\_|_| |_|\__, |
##                               |___/                                                                        |___/ 
##   
## Match them to the collection of "tkTracks" which are defined in the common cff
staToTkMatch = cms.EDProducer("MatcherUsingTracks",
    src     = cms.InputTag("staTracks"), # all standalone muons
    matched = cms.InputTag("tkTracks"),  # to all tk tracks
    algorithm = cms.string("byDirectComparison"), # using parameters at PCA
    srcTrack = cms.string("tracker"),  # 'staTracks' is a 'RecoChargedCandidate', so it thinks
    srcState = cms.string("atVertex"), # it has a 'tracker' track, not a standalone one
    matchedTrack = cms.string("tracker"),
    matchedState = cms.string("atVertex"),
    maxDeltaR        = cms.double(1.),   # large range in DR
    maxDeltaEta      = cms.double(0.4),  # small in eta, which is more precise
    maxDeltaLocalPos = cms.double(100),
    maxDeltaPtRel    = cms.double(5),
    sortBy           = cms.string("deltaR"),
)
staPassingTk = cms.EDProducer("MatchedCandidateSelector",
    src   = cms.InputTag("staProbes"),
    match = cms.InputTag("staToTkMatch"),
)
staToTkMatchNoJPsi     = staToTkMatch.clone(matched = 'tkTracksNoJPsi')
staToTkMatchNoBestJPsi = staToTkMatch.clone(matched = 'tkTracksNoBestJPsi')
staPassingTkNoJPsi     = staPassingTk.clone(match = 'staToTkMatchNoJPsi')
staPassingTkNoBestJPsi = staPassingTk.clone(match = 'staToTkMatchNoBestJPsi')

allProbesTracking = cms.Sequence(
    staTracks * staProbes 
)
allPassingProbesTracking = cms.Sequence(
    staToTkMatch           * staPassingTk           +
    staToTkMatchNoJPsi     * staPassingTkNoJPsi     +
    staToTkMatchNoBestJPsi * staPassingTkNoBestJPsi  
)

##    _____ ___   ____    ____       _          
##   |_   _( _ ) |  _ \  |  _ \ __ _(_)_ __ ___ 
##     | | / _ \/\ |_) | | |_) / _` | | '__/ __|
##     | || (_>  <  __/  |  __/ (_| | | |  \__ \
##     |_| \___/\/_|     |_|   \__,_|_|_|  |___/
##                                              
##   
tpGlbSta = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("tagMuons2Mu@+ staProbes@-"), # charge coniugate states are implied
    cut   = cms.string("%f < mass < %f" % MASS_RANGE_STA),
)

allTPPairsTracking = cms.Sequence(
    tpGlbSta
)


staMcMatch = muMcMatch.clone(src = "staTracks", distMin = 0.6)
allMcMatchesTracking = cms.Sequence(staMcMatch)


##    _____           _       _ ____            _            _   _ _____            _      
##   |_   _|_ _  __ _( )_ __ ( )  _ \ _ __ ___ | |__   ___  | \ | |_   _|   _ _ __ | | ___ 
##     | |/ _` |/ _` |/| '_ \|/| |_) | '__/ _ \| '_ \ / _ \ |  \| | | || | | | '_ \| |/ _ \
##     | | (_| | (_| | | | | | |  __/| | | (_) | |_) |  __/ | |\  | | || |_| | |_) | |  __/
##     |_|\__,_|\__, | |_| |_| |_|   |_|  \___/|_.__/ \___| |_| \_| |_| \__,_| .__/|_|\___|
##              |___/                                                        |_|           
##
histoTracking = tnpTreeProducer.clone(
    tagProbePairs = cms.InputTag("tpGlbSta"),
    flags = cms.PSet(
        # Passing definition
        passing = cms.InputTag('staPassingTk'),
        passingNoJPsi     = cms.InputTag('staPassingTkNoJPsi'),
        passingNoBestJPsi = cms.InputTag('staPassingTkNoBestJPsi'),
        # Quality cuts
        hasValidHits = cms.string('track.numberOfValidHits > 0'),
    ),
    # MC Matching configurables
    probeMatches  = cms.InputTag("staMcMatch"),
    allProbes = cms.InputTag("staProbes"),
)
histoTracking.variables.match_deltaR   = cms.InputTag("staToTkMatch", "deltaR")
histoTracking.variables.match_deltaEta = cms.InputTag("staToTkMatch", "deltaEta")
histoTracking.variables.match_deltaR_NoJPsi   = cms.InputTag("staToTkMatchNoJPsi", "deltaR")
histoTracking.variables.match_deltaEta_NoJPsi = cms.InputTag("staToTkMatchNoJPsi", "deltaEta")
histoTracking.variables.match_deltaR_NoBestJPsi   = cms.InputTag("staToTkMatchNoBestJPsi", "deltaR")
histoTracking.variables.match_deltaEta_NoBestJPsi = cms.InputTag("staToTkMatchNoBestJPsi", "deltaEta")

allTPHistosTracking = cms.Sequence(
        histoTracking     
)

tnpSequenceTracking = cms.Sequence(
    allProbesTracking *
    (tkTracksNoJPsi + tkTracksNoBestJPsi) *
    allPassingProbesTracking *
    allTPPairsTracking   *
    allMcMatchesTracking * 
    allTPHistosTracking
)

import FWCore.ParameterSet.Config as cms

from MuonAnalysis.TagAndProbe.jpsi.tp_from_skim_common_cff import *

##    ____  _                  _    _    _                    ____            _               
##   / ___|| |_ __ _ _ __   __| |  / \  | | ___  _ __   ___  |  _ \ _ __ ___ | |__   ___  ___ 
##   \___ \| __/ _` | '_ \ / _` | / _ \ | |/ _ \| '_ \ / _ \ | |_) | '__/ _ \| '_ \ / _ \/ __|
##    ___) | || (_| | | | | (_| |/ ___ \| | (_) | | | |  __/ |  __/| | | (_) | |_) |  __/\__ \
##   |____/ \__\__,_|_| |_|\__,_/_/   \_\_|\___/|_| |_|\___| |_|   |_|  \___/|_.__/ \___||___/
##                                                                                            
##   
staProbes = cms.EDFilter("CandViewRefSelector",
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
justhpTkTracksNoJPsi     = tkTracksNoJPsi.clone(src = 'justhpTkTracks')
justhpTkTracksNoBestJPsi = tkTracksNoBestJPsi.clone(src = 'justhpTkTracks')

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

staToHpTkMatch           = staToTkMatch.clone(matched = 'justhpTkTracks')
staToHpTkMatchNoJPsi     = staToTkMatch.clone(matched = 'justhpTkTracksNoJPsi')
staToHpTkMatchNoBestJPsi = staToTkMatch.clone(matched = 'justhpTkTracksNoBestJPsi')
staPassingHpTk           = staPassingTk.clone(match = 'staToHpTkMatch')
staPassingHpTkNoJPsi     = staPassingTk.clone(match = 'staToHpTkMatchNoJPsi')
staPassingHpTkNoBestJPsi = staPassingTk.clone(match = 'staToHpTkMatchNoBestJPsi')

allProbesTracking = cms.Sequence(
    staTracks * staProbes 
)
allPassingProbesTracking = cms.Sequence(
    staToHpTkMatch           * staPassingHpTk           +
    staToHpTkMatchNoJPsi     * staPassingHpTkNoJPsi     +
    staToHpTkMatchNoBestJPsi * staPassingHpTkNoBestJPsi +
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

##    _____      _               ___        __       
##   | ____|_  _| |_ _ __ __ _  |_ _|_ __  / _| ___  
##   |  _| \ \/ / __| '__/ _` |  | || '_ \| |_ / _ \ 
##   | |___ >  <| |_| | | (_| |  | || | | |  _| (_) |
##   |_____/_/\_\\__|_|  \__,_| |___|_| |_|_|  \___/ 
##                                                   
##   
countNearbyTracks4Tracking = cms.EDProducer("NearbyCandCountComputer",
    probes  = cms.InputTag("staTracks"),
    objects = cms.InputTag("justhpTkTracks"),
    deltaR  = cms.double(1.0)
)

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
histoTracking.variables.ntracksDR1 = cms.InputTag("countNearbyTracks4Tracking")

histoTrackingHp = histoTracking.clone()
histoTrackingHp.flags.passing           = cms.InputTag('staPassingTk')
histoTrackingHp.flags.passingNoJPsi     = cms.InputTag('staPassingTkNoJPsi')
histoTrackingHp.flags.passingNoBestJPsi = cms.InputTag('staPassingTkNoBestJPsi')
histoTrackingHp.variables.match_deltaR   = cms.InputTag("staToHpTkMatch", "deltaR")
histoTrackingHp.variables.match_deltaEta = cms.InputTag("staToHpTkMatch", "deltaEta")
histoTrackingHp.variables.match_deltaR_NoJPsi   = cms.InputTag("staToHpTkMatchNoJPsi", "deltaR")
histoTrackingHp.variables.match_deltaEta_NoJPsi = cms.InputTag("staToHpTkMatchNoJPsi", "deltaEta")
histoTrackingHp.variables.match_deltaR_NoBestJPsi   = cms.InputTag("staToHpTkMatchNoBestJPsi", "deltaR")
histoTrackingHp.variables.match_deltaEta_NoBestJPsi = cms.InputTag("staToHpTkMatchNoBestJPsi", "deltaEta")
histoTrackingHp.variables.ntracksDR1 = cms.InputTag("countNearbyTracks4Tracking")

allTPHistosTracking = cms.Sequence(
        histoTracking     
        + histoTrackingHp     
)

##    _____               _       ___              _ _ _         
##   |_   _| __ __ _  ___| | __  / _ \ _   _  __ _| (_) |_ _   _ 
##     | || '__/ _` |/ __| |/ / | | | | | | |/ _` | | | __| | | |
##     | || | | (_| | (__|   <  | |_| | |_| | (_| | | | |_| |_| |
##     |_||_|  \__,_|\___|_|\_\  \__\_\\__,_|\__,_|_|_|\__|\__, |
##                                                         |___/ 
##    _____  __  __ _      _                 _           
##   | ____|/ _|/ _(_) ___(_) ___ _ __   ___(_) ___  ___ 
##   |  _| | |_| |_| |/ __| |/ _ \ '_ \ / __| |/ _ \/ __|
##   | |___|  _|  _| | (__| |  __/ | | | (__| |  __/\__ \
##   |_____|_| |_| |_|\___|_|\___|_| |_|\___|_|\___||___/
##                                                       
##   

muonsNoTrackQualityCuts = cms.EDFilter("PATMuonRefSelector",
    src = cms.InputTag("patMuons"),
    cut = cms.string("isGlobalMuon || muonID('TMLastStationAngTight')"),
)
tpGlbMuNoTkQ = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("tagMuons2Mu@+ muonsNoTrackQualityCuts@-"), # charge coniugate states are implied
    cut   = cms.string("%f < mass < %f" % MASS_RANGE_STA),
)
histoTrackQuality = tnpTreeProducer.clone(
    tagProbePairs = cms.InputTag("tpGlbMuNoTkQ"),
    flags = cms.PSet(
        # Passing definition
        highPurity  = cms.string("track.quality('highPurity')"),
        QTF_Tracks  = cms.string(TRACK_CUTS),
        VBTF_Tracks = cms.string("track.numberOfValidHits > 10 && track.hitPattern.numberOfValidPixelHits > 0 && abs(dB) < 0.2"),
        Hits12      = cms.string("track.numberOfValidHits >= 12"),
        PxlLay2     = cms.string("track.hitPattern.pixelLayersWithMeasurement >= 2"),
        Chi2Ndf4    = cms.string("track.normalizedChi2 < 4"),
    ),
    # MC Matching configurables
    probeMatches  = cms.InputTag("muMcMatch"),
    allProbes = cms.InputTag("muonsNoTrackQualityCuts"),
)
histoTrackQuality.variables.hits = cms.string("track.numberOfValidHits")
histoTrackQuality.variables.pixelHits = cms.string("track.hitPattern.numberOfValidPixelHits")
histoTrackQuality.variables.pixelLayers = cms.string("track.hitPattern.pixelLayersWithMeasurement")
histoTrackQuality.variables.expHitsOut = cms.string("track.trackerExpectedHitsOuter.numberOfLostHits")
histoTrackQuality.variables.expHitsIn  = cms.string("track.trackerExpectedHitsInner.numberOfLostHits")
histoTrackQuality.variables.dB    = cms.string("dB")
histoTrackQuality.variables.chi2n = cms.string("track.normalizedChi2")

tnpSequenceTrackQuality = cms.Sequence(
    muonsNoTrackQualityCuts +
    tpGlbMuNoTkQ +
    histoTrackQuality
)

tnpSequenceTracking = cms.Sequence(
    allProbesTracking *
    (tkTracksNoJPsi + tkTracksNoBestJPsi +
     justhpTkTracksNoJPsi + justhpTkTracksNoBestJPsi) *
    allPassingProbesTracking *
    allTPPairsTracking   *
    allMcMatchesTracking * 
    countNearbyTracks4Tracking *
    allTPHistosTracking +
    tnpSequenceTrackQuality
)

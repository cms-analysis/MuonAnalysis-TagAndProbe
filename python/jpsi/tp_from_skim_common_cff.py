import FWCore.ParameterSet.Config as cms

##   __     __         _       _     _                             _    ____      _       
##   \ \   / /_ _ _ __(_) __ _| |__ | | ___  ___    __ _ _ __   __| |  / ___|   _| |_ ___ 
##    \ \ / / _` | '__| |/ _` | '_ \| |/ _ \/ __|  / _` | '_ \ / _` | | |  | | | | __/ __|
##     \ V / (_| | |  | | (_| | |_) | |  __/\__ \ | (_| | | | | (_| | | |__| |_| | |_\__ \
##      \_/ \__,_|_|  |_|\__,_|_.__/|_|\___||___/  \__,_|_| |_|\__,_|  \____\__,_|\__|___/
##                                                                                        
##   

MASS_RANGE = ( 2.8, 3.5 )
MASS_RANGE_STA = ( 2, 5 )
TRACK_CUTS = ("track.numberOfValidHits > 11 && track.hitPattern.pixelLayersWithMeasurement > 1 "+
              "&& track.normalizedChi2 < 5 "+
              "&& abs(track.d0) < 2 && abs(track.dz) < 30")
PT_ETA_CUTS = "(pt > 3 || (abs(eta)>1 && p > 2.6))" ## the enclosing () are very important, because there's an "||"
PASSING_GLB_CUT = "isGlobalMuon && globalTrack.normalizedChi2 < 20"
PASSING_TM_CUT = "!("+PASSING_GLB_CUT+") && isTrackerMuon && muonID('TMLastStationAngTight')";

PASS_HLT_1MU = "!triggerObjectMatchesByFilter('hltSingleMu3L3Filtered3').empty()"
PASS_HLT_2MU = "!triggerObjectMatchesByFilter('hltDoubleMuLevel1PathL1OpenFiltered').empty()"
TAG_CUTS_1MU = PASSING_GLB_CUT + " && " + TRACK_CUTS +' && '+ PASS_HLT_1MU
TAG_CUTS_2MU = PASSING_GLB_CUT + " && " + TRACK_CUTS +' && '+ PASS_HLT_2MU

##     ____ _       _           _ __  __                     _____               
##    / ___| | ___ | |__   __ _| |  \/  |_   _  ___  _ __   |_   _|_ _  __ _ ___ 
##   | |  _| |/ _ \| '_ \ / _` | | |\/| | | | |/ _ \| '_ \    | |/ _` |/ _` / __|
##   | |_| | | (_) | |_) | (_| | | |  | | |_| | (_) | | | |   | | (_| | (_| \__ \
##    \____|_|\___/|_.__/ \__,_|_|_|  |_|\__,_|\___/|_| |_|   |_|\__,_|\__, |___/
##                                                                     |___/     
##   
tagMuons1Mu = cms.EDFilter("PATMuonRefSelector",
    src = cms.InputTag("patMuons"),
    cut = cms.string(TAG_CUTS_1MU), 
)
tagMuons2Mu = cms.EDFilter("PATMuonRefSelector",
    src = cms.InputTag("patMuons"),
    cut = cms.string(TAG_CUTS_2MU), 
)

##     ____                 _   _                  _        
##    / ___| ___   ___   __| | | |_ _ __ __ _  ___| | _____ 
##   | |  _ / _ \ / _ \ / _` | | __| '__/ _` |/ __| |/ / __|
##   | |_| | (_) | (_) | (_| | | |_| | | (_| | (__|   <\__ \
##    \____|\___/ \___/ \__,_|  \__|_|  \__,_|\___|_|\_\___/
##                                                          
##   
betterTracks = cms.EDFilter("TrackSelector",
    src = cms.InputTag("goodTracks"),
    cut = cms.string(TRACK_CUTS.replace("track.","")), # this is a Track object, so I have to remove the 'track.'
)
tkTracks  = cms.EDProducer("ConcreteChargedCandidateProducer", 
    src  = cms.InputTag("betterTracks"),      
    particleType = cms.string("mu+"),
) 

##    __  __  ____   __  __       _       _               
##   |  \/  |/ ___| |  \/  | __ _| |_ ___| |__   ___  ___ 
##   | |\/| | |     | |\/| |/ _` | __/ __| '_ \ / _ \/ __|
##   | |  | | |___  | |  | | (_| | || (__| | | |  __/\__ \
##   |_|  |_|\____| |_|  |_|\__,_|\__\___|_| |_|\___||___/
##                                                        
muMcMatch = cms.EDFilter("MCTruthDeltaRMatcherNew",
    pdgId = cms.vint32(13),
    src = cms.InputTag("patMuons"),
    distMin = cms.double(0.3),
    matched = cms.InputTag("genMuons")
)

##    _____           _       _ ____            _            _   _ _____            _      
##   |_   _|_ _  __ _( )_ __ ( )  _ \ _ __ ___ | |__   ___  | \ | |_   _|   _ _ __ | | ___ 
##     | |/ _` |/ _` |/| '_ \|/| |_) | '__/ _ \| '_ \ / _ \ |  \| | | || | | | '_ \| |/ _ \
##     | | (_| | (_| | | | | | |  __/| | | (_) | |_) |  __/ | |\  | | || |_| | |_) | |  __/
##     |_|\__,_|\__, | |_| |_| |_|   |_|  \___/|_.__/ \___| |_| \_| |_| \__,_| .__/|_|\___|
##              |___/                                                        |_|           
##

tnpTreeProducer = cms.EDAnalyzer("TagProbeFitTreeProducer",
    # choice of tag and probe pairs, and arbitration
    tagProbePairs = cms.InputTag("REPLACE_ME"),
    arbitration   = cms.string("OneProbe"),
    # probe variables 
    variables = cms.PSet(
        pt  = cms.string("pt"),
        p   = cms.string("p"),
        eta = cms.string("eta"),
        phi = cms.string("phi"),
    ),
    # choice of what defines a 'passing' probe
    flags = cms.PSet(),
    ## variables relative to the _tag_
    tagVariables = cms.PSet(
        pt  = cms.string("pt"),
        eta = cms.string("eta"),
     ),
    tagFlags = cms.PSet(
        HLTMu3         = cms.string("!triggerObjectMatchesByFilter('hltSingleMu3L3Filtered3').empty()"),
        L1SingleMuOpen = cms.string("!triggerObjectMatchesByFilter('hltL1MuOpenL1Filtered0').empty()"),
    ),
    ## MC-related info
    isMC = cms.bool(True),
    makeMCUnbiasTree = cms.bool(False),       ## NO! 'unbias' efficiency on a skim is
    checkMotherInUnbiasEff = cms.bool(True),  ##      biased _a lot_ by the skim tags
    tagMatches = cms.InputTag("muMcMatch"),
    motherPdgId = cms.int32(443),
    ## These two MC things depend on the specific choice of probes
    probeMatches  = cms.InputTag("REPLACE_ME"),
    allProbes     = cms.InputTag("REPLACE_ME"), 
)


tnpCommonSequence = cms.Sequence(
    betterTracks * tkTracks +
    tagMuons1Mu + 
    tagMuons2Mu + 
    muMcMatch 
)


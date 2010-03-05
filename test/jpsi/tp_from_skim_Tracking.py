import FWCore.ParameterSet.Config as cms

process = cms.Process("TagProbe")

process.load('FWCore.MessageService.MessageLogger_cfi')
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.MessageLogger.cerr.FwkReport.reportEvery = 100

massRange = ( 2.0, 5.0 )
TRACK_CUTS="track.numberOfValidHits >= 10 && track.normalizedChi2 < 5 && abs(track.d0) < 2 && abs(track.dz) < 30"
PT_ETA_CUTS="p > 2" # currently, it looks like there is no reason to put any tighter cut here 


process.source = cms.Source("PoolSource", 
    fileNames = cms.untracked.vstring(
        'root://pcmssd12.cern.ch//data/gpetrucc/tnp/ppMuX_skimJPsi_1.root',
        'root://pcmssd12.cern.ch//data/gpetrucc/tnp/ppMuX_skimJPsi_2.root',
        'root://pcmssd12.cern.ch//data/gpetrucc/tnp/JPsiMuMu_skimJPsi_1.root',
    )
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )    

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Geometry.CommonDetUnit.globalTrackingGeometry_cfi")
process.load("TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorAny_cfi")
process.load("TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorAlong_cfi")
process.load("TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorOpposite_cfi")
process.GlobalTag.globaltag = cms.string('MC_3XY_V18::All')

##    ____  _                  _    _    _                    ____            _               
##   / ___|| |_ __ _ _ __   __| |  / \  | | ___  _ __   ___  |  _ \ _ __ ___ | |__   ___  ___ 
##   \___ \| __/ _` | '_ \ / _` | / _ \ | |/ _ \| '_ \ / _ \ | |_) | '__/ _ \| '_ \ / _ \/ __|
##    ___) | || (_| | | | | (_| |/ ___ \| | (_) | | | |  __/ |  __/| | | (_) | |_) |  __/\__ \
##   |____/ \__\__,_|_| |_|\__,_/_/   \_\_|\___/|_| |_|\___| |_|   |_|  \___/|_.__/ \___||___/
##                                                                                            
##   
process.staTracks = cms.EDProducer("ConcreteChargedCandidateProducer", 
    src  = cms.InputTag("standAloneMuons","UpdatedAtVtx"), 
    particleType = cms.string("mu+"),
)
process.staProbes = cms.EDProducer("CandViewRefSelector",
    src = cms.InputTag("staTracks"),
    cut = cms.string(PT_ETA_CUTS),
)

##    _____               
##   |_   _|_ _  __ _ ___ 
##     | |/ _` |/ _` / __|
##     | | (_| | (_| \__ \
##     |_|\__,_|\__, |___/
##              |___/     
##   
PASS_HLT = "!triggerObjectMatchesByFilter('hltDoubleMuLevel1PathL1OpenFiltered').empty()";
process.tagMuons = cms.EDFilter("PATMuonRefSelector",
    src = cms.InputTag("patMuons"),
    cut = cms.string("isGlobalMuon && " + TRACK_CUTS + " && " +PASS_HLT ), 
    filter = cms.bool(True),
)

##    _____               _           __              __  __       _       _     _             
##   |_   _| __ __ _  ___| | _____   / _| ___  _ __  |  \/  | __ _| |_ ___| |__ (_)_ __   __ _ 
##     | || '__/ _` |/ __| |/ / __| | |_ / _ \| '__| | |\/| |/ _` | __/ __| '_ \| | '_ \ / _` |
##     | || | | (_| | (__|   <\__ \ |  _| (_) | |    | |  | | (_| | || (__| | | | | | | | (_| |
##     |_||_|  \__,_|\___|_|\_\___/ |_|  \___/|_|    |_|  |_|\__,_|\__\___|_| |_|_|_| |_|\__, |
##                                                                                       |___/ 
##   
## Enforce the quality cuts, because we consider them "tracking" in our factorization of efficiencies
process.betterTracks = cms.EDFilter("TrackSelector",
    src = cms.InputTag("goodTracks"),
    cut = cms.string(TRACK_CUTS.replace("track.","")), # this is a Track object, so I have to remove the 'track.'
)
process.tkTracks  = cms.EDProducer("ConcreteChargedCandidateProducer", 
    src  = cms.InputTag("betterTracks"),      
    particleType = cms.string("mu+"),
) 



##    ____               _               ____            _                   _____               _    _             
##   |  _ \ __ _ ___ ___(_)_ __   __ _  |  _ \ _ __ ___ | |__   ___  ___ _  |_   _| __ __ _  ___| | _(_)_ __   __ _ 
##   | |_) / _` / __/ __| | '_ \ / _` | | |_) | '__/ _ \| '_ \ / _ \/ __(_)   | || '__/ _` |/ __| |/ / | '_ \ / _` |
##   |  __/ (_| \__ \__ \ | | | | (_| | |  __/| | | (_) | |_) |  __/\__ \_    | || | | (_| | (__|   <| | | | | (_| |
##   |_|   \__,_|___/___/_|_| |_|\__, | |_|   |_|  \___/|_.__/ \___||___(_)   |_||_|  \__,_|\___|_|\_\_|_| |_|\__, |
##                               |___/                                                                        |___/ 
##   
process.staToTkMatch = cms.EDProducer("MatcherUsingTracks",
    src     = cms.InputTag("staTracks"), # all standalone muons
    matched = cms.InputTag("tkTracks"),  # to all tk tracks
    algorithm = cms.string("byDirectComparison"), # using parameters at PCA
    srcTrack = cms.string("tracker"),  # 'staTracks' is a 'RecoChargedCandidate', so it thinks
    srcState = cms.string("atVertex"), # it has a 'tracker' track, not a standalone one
    matchedTrack = cms.string("tracker"),
    matchedState = cms.string("atVertex"),
    maxDeltaR        = cms.double(1.),   # large range in DR
    maxDeltaEta      = cms.double(0.3),  # small in eta, which is more precise
    maxDeltaLocalPos = cms.double(100),
    maxDeltaPtRel    = cms.double(3),
    sortBy           = cms.string("deltaR"),
)
process.staPassingTk = cms.EDProducer("MatchedCandidateSelector",
    src   = cms.InputTag("staProbes"),
    match = cms.InputTag("staToTkMatch"),
)

process.allTagsAndProbes = cms.Sequence(
    process.tagMuons +
    process.betterTracks * process.tkTracks +
    process.staTracks    * process.staProbes +
    process.staToTkMatch * process.staPassingTk 
)

##    _____ ___   ____    ____       _          
##   |_   _( _ ) |  _ \  |  _ \ __ _(_)_ __ ___ 
##     | | / _ \/\ |_) | | |_) / _` | | '__/ __|
##     | || (_>  <  __/  |  __/ (_| | | |  \__ \
##     |_| \___/\/_|     |_|   \__,_|_|_|  |___/
##                                              
##   
process.tpGlbSta = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("tagMuons@+ staProbes@-"), # charge coniugate states are implied
    cut   = cms.string("%f < mass < %f" % massRange),
)

##    __  __  ____   __  __       _       _               
##   |  \/  |/ ___| |  \/  | __ _| |_ ___| |__   ___  ___ 
##   | |\/| | |     | |\/| |/ _` | __/ __| '_ \ / _ \/ __|
##   | |  | | |___  | |  | | (_| | || (__| | | |  __/\__ \
##   |_|  |_|\____| |_|  |_|\__,_|\__\___|_| |_|\___||___/
##                                                        
process.tagMcMatch = cms.EDFilter("MCTruthDeltaRMatcherNew",
    pdgId = cms.vint32(13),
    src = cms.InputTag("patMuons"),
    distMin = cms.double(0.3),
    matched = cms.InputTag("genMuons")
)
process.staMcMatch = process.tagMcMatch.clone(src = "staTracks", distMin = 0.6)

process.allMcMatches = cms.Sequence(process.tagMcMatch + process.staMcMatch)


##    _____           _       _ ____            _            _   _ _____            _      
##   |_   _|_ _  __ _( )_ __ ( )  _ \ _ __ ___ | |__   ___  | \ | |_   _|   _ _ __ | | ___ 
##     | |/ _` |/ _` |/| '_ \|/| |_) | '__/ _ \| '_ \ / _ \ |  \| | | || | | | '_ \| |/ _ \
##     | | (_| | (_| | | | | | |  __/| | | (_) | |_) |  __/ | |\  | | || |_| | |_) | |  __/
##     |_|\__,_|\__, | |_| |_| |_|   |_|  \___/|_.__/ \___| |_| \_| |_| \__,_| .__/|_|\___|
##              |___/                                                        |_|           
##
process.histo = cms.EDAnalyzer("TagProbeFitTreeProducer",
    variables = cms.PSet(
        pt  = cms.string("pt"),
        p   = cms.string("p"),
        eta = cms.string("eta"),
        phi = cms.string("phi"),
    ),
    # choice of tag and probe pairs, and arbitration
    tagProbePairs = cms.InputTag("tpGlbSta"),
    arbitration   = cms.string("OneProbe"),
    # choice of what defines a 'passing' probe
    flags = cms.PSet(passing = cms.InputTag('staPassingTk')),
    isMC = cms.bool(True),
    tagMatches = cms.InputTag("tagMcMatch"),
    probeMatches  = cms.InputTag("staMcMatch"),
    motherPdgId = cms.int32(443),
    makeMCUnbiasTree = cms.bool(False),        ## NO! 'unbias' efficiency on a skim is
    #checkMotherInUnbiasEff = cms.bool(True),  ##      biased _a lot_ by the skim tags
    #allProbes = cms.InputTag("staProbes"),
)



##    ____       _   _     
##   |  _ \ __ _| |_| |__  
##   | |_) / _` | __| '_ \ 
##   |  __/ (_| | |_| | | |
##   |_|   \__,_|\__|_| |_|
##                         
process.tagAndProbe = cms.Path( 
    process.allTagsAndProbes *
    process.tpGlbSta   *
    process.allMcMatches * 
    process.histo
)

##     ___        _               _     _   _ _     _            
##    / _ \ _   _| |_ _ __  _   _| |_  | | | (_)___| |_ ___  ___ 
##   | | | | | | | __| '_ \| | | | __| | |_| | / __| __/ _ \/ __|
##   | |_| | |_| | |_| |_) | |_| | |_  |  _  | \__ \ || (_) \__ \
##    \___/ \__,_|\__| .__/ \__,_|\__| |_| |_|_|___/\__\___/|___/
##                   |_|                                         
##   
process.TFileService = cms.Service("TFileService", fileName = cms.string("tnpJPsi_Tracking.root"))



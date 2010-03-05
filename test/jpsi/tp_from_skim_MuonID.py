import FWCore.ParameterSet.Config as cms

process = cms.Process("TagProbe")

process.load('FWCore.MessageService.MessageLogger_cfi')
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.MessageLogger.cerr.FwkReport.reportEvery = 100

##   __     __         _       _     _                             _    ____      _       
##   \ \   / /_ _ _ __(_) __ _| |__ | | ___  ___    __ _ _ __   __| |  / ___|   _| |_ ___ 
##    \ \ / / _` | '__| |/ _` | '_ \| |/ _ \/ __|  / _` | '_ \ / _` | | |  | | | | __/ __|
##     \ V / (_| | |  | | (_| | |_) | |  __/\__ \ | (_| | | | | (_| | | |__| |_| | |_\__ \
##      \_/ \__,_|_|  |_|\__,_|_.__/|_|\___||___/  \__,_|_| |_|\__,_|  \____\__,_|\__|___/
##                                                                                        
##   
MASS_RANGE = ( 2.8, 3.5 )
TRACK_CUTS = ("track.numberOfValidHits > 11 && track.hitPattern.pixelLayersWithMeasurement > 1 "+
              "&& track.normalizedChi2 < 5 "+
              "&& abs(track.d0) < 2 && abs(track.dz) < 30")
PT_ETA_CUTS = "(pt > 3 || (abs(eta)>1 && p > 2.6))" ## the enclosing () are very important, because there's an "||"
PASSING_GLB_CUT = "isGlobalMuon && globalTrack.normalizedChi2 < 20"
PASSING_TM_CUT = "!("+PASSING_GLB_CUT+") && isTrackerMuon && muonID('TMLastStationAngTight')";

PASS_HLT = "!triggerObjectMatchesByFilter('hltSingleMu3L3Filtered3').empty()"
TAG_CUTS = PASSING_GLB_CUT + " && " + TRACK_CUTS +' && '+ PASS_HLT 

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

##    ____                   _____               _      ____            _               
##   | __ )  __ _ _ __ ___  |_   _| __ __ _  ___| | __ |  _ \ _ __ ___ | |__   ___  ___ 
##   |  _ \ / _` | '__/ _ \   | || '__/ _` |/ __| |/ / | |_) | '__/ _ \| '_ \ / _ \/ __|
##   | |_) | (_| | | |  __/   | || | | (_| | (__|   <  |  __/| | | (_) | |_) |  __/\__ \
##   |____/ \__,_|_|  \___|   |_||_|  \__,_|\___|_|\_\ |_|   |_|  \___/|_.__/ \___||___/
##                                                                                      
##   
process.betterTracks = cms.EDFilter("TrackSelector",
    src = cms.InputTag("goodTracks"),
    cut = cms.string(TRACK_CUTS.replace("track.","")), # this is a Track object, so I have to remove the 'track.'
)
process.tkTracks  = cms.EDProducer("ConcreteChargedCandidateProducer", 
    src  = cms.InputTag("betterTracks"),      
    particleType = cms.string("mu+"),
) 
process.tkProbes = cms.EDProducer("CandViewRefSelector",
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
process.calProbes = cms.EDFilter("PATMuonRefSelector",
    src = cms.InputTag("patMuons"),
    cut = cms.string("track.isNonnull && isCaloMuon && "+PT_ETA_CUTS+" && "+TRACK_CUTS),
)

process.allProbes = cms.Sequence(
    process.betterTracks * process.tkTracks * process.tkProbes +
    process.calProbes 
)

##     ____ _       _           _ __  __                     _____               
##    / ___| | ___ | |__   __ _| |  \/  |_   _  ___  _ __   |_   _|_ _  __ _ ___ 
##   | |  _| |/ _ \| '_ \ / _` | | |\/| | | | |/ _ \| '_ \    | |/ _` |/ _` / __|
##   | |_| | | (_) | |_) | (_| | | |  | | |_| | (_) | | | |   | | (_| | (_| \__ \
##    \____|_|\___/|_.__/ \__,_|_|_|  |_|\__,_|\___/|_| |_|   |_|\__,_|\__, |___/
##                                                                     |___/     
##   
process.tagMuons = cms.EDFilter("PATMuonRefSelector",
    src = cms.InputTag("patMuons"),
    cut = cms.string(TAG_CUTS), 
    filter = cms.bool(True),
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
process.muonsGlb = cms.EDFilter("PATMuonRefSelector",
    src = cms.InputTag("patMuons"),
    cut = cms.string(PASSING_GLB_CUT),
)
process.muonsTM = process.muonsGlb.clone(cut = PASSING_TM_CUT)

process.tkToGlbMatch = cms.EDProducer("MatcherUsingTracks",
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
process.tkPassingGlb = cms.EDProducer("MatchedCandidateSelector",
    src   = cms.InputTag("tkProbes"),
    match = cms.InputTag("tkToGlbMatch"),
)
process.tkToTMMatch = process.tkToGlbMatch.clone(matched = 'muonsTM')
process.tkPassingTM = process.tkPassingGlb.clone( match = 'tkToTMMatch')

process.tkPassingProbes = cms.Sequence(
    process.muonsGlb * process.tkToGlbMatch * process.tkPassingGlb +
    process.muonsTM  * process.tkToTMMatch  * process.tkPassingTM 
)
##    _____ ___   ____    ____       _          
##   |_   _( _ ) |  _ \  |  _ \ __ _(_)_ __ ___ 
##     | | / _ \/\ |_) | | |_) / _` | | '__/ __|
##     | || (_>  <  __/  |  __/ (_| | | |  \__ \
##     |_| \___/\/_|     |_|   \__,_|_|_|  |___/
##                                              
##   
process.tpGlbTk = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("tagMuons@+ tkProbes@-"), # charge coniugate states are implied
    cut   = cms.string("%f < mass < %f" % MASS_RANGE),
)
process.tpGlbCal = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("tagMuons@+ calProbes@-"), # charge coniugate states are implied
    cut   = cms.string("%f < mass < %f" % MASS_RANGE),
)
process.allTPPairs = cms.Sequence(process.tpGlbTk + process.tpGlbCal)


##    __  __  ____   __  __       _       _               
##   |  \/  |/ ___| |  \/  | __ _| |_ ___| |__   ___  ___ 
##   | |\/| | |     | |\/| |/ _` | __/ __| '_ \ / _ \/ __|
##   | |  | | |___  | |  | | (_| | || (__| | | |  __/\__ \
##   |_|  |_|\____| |_|  |_|\__,_|\__\___|_| |_|\___||___/
##                                                        
process.muMcMatch = cms.EDFilter("MCTruthDeltaRMatcherNew",
    pdgId = cms.vint32(13),
    src = cms.InputTag("patMuons"),
    distMin = cms.double(0.3),
    matched = cms.InputTag("genMuons")
)
process.tkMcMatch  = process.muMcMatch.clone(src = "tkTracks")

process.allMcMatches = cms.Sequence(process.muMcMatch + process.tkMcMatch)


##    _____           _       _ ____            _            _   _ _____            _      
##   |_   _|_ _  __ _( )_ __ ( )  _ \ _ __ ___ | |__   ___  | \ | |_   _|   _ _ __ | | ___ 
##     | |/ _` |/ _` |/| '_ \|/| |_) | '__/ _ \| '_ \ / _ \ |  \| | | || | | | '_ \| |/ _ \
##     | | (_| | (_| | | | | | |  __/| | | (_) | |_) |  __/ | |\  | | || |_| | |_) | |  __/
##     |_|\__,_|\__, | |_| |_| |_|   |_|  \___/|_.__/ \___| |_| \_| |_| \__,_| .__/|_|\___|
##              |___/                                                        |_|           
##
recoCommonStuff = cms.PSet(
    variables = cms.PSet(
        pt  = cms.string("pt"),
        p   = cms.string("p"),
        eta = cms.string("eta"),
        phi = cms.string("phi"),
    )
)
mcTruthCommonStuff = cms.PSet(
    isMC = cms.bool(True),
    makeMCUnbiasTree = cms.bool(False),        ## NO! 'unbias' efficiency on a skim is
    #checkMotherInUnbiasEff = cms.bool(True),  ##      biased _a lot_ by the skim tags
    tagMatches = cms.InputTag("muMcMatch"),
    motherPdgId = cms.int32(443),
)

#####
## Mu from Tk
process.histoMuFromTk = cms.EDAnalyzer("TagProbeFitTreeProducer",
    ## pick the defaults
    recoCommonStuff, mcTruthCommonStuff,
    # choice of tag and probe pairs, and arbitration
    tagProbePairs = cms.InputTag("tpGlbTk"),
    arbitration   = cms.string("OneProbe"),
    # choice of what defines a 'passing' probe
    flags = cms.PSet(
        Glb = cms.InputTag("tkPassingGlb"),
        TM  = cms.InputTag("tkPassingTM"),
    ),
    ## These two MC things depend on the specific choice of probes
    probeMatches  = cms.InputTag("tkMcMatch"),
    #allProbes     = cms.InputTag("tkProbes"), # NO 'unbias' efficiency on skims
)

#####
## Mu from Cal
process.histoMuFromCal = cms.EDAnalyzer("TagProbeFitTreeProducer",
    ## pick the defaults
    recoCommonStuff, mcTruthCommonStuff,
    # choice of tag and probe pairs, and arbitration
    tagProbePairs = cms.InputTag("tpGlbCal"),
    arbitration   = cms.string("OneProbe"),
    # choice of what defines a 'passing' probe
    flags = cms.PSet(
        Glb = cms.string(PASSING_GLB_CUT),
        TM  = cms.string(PASSING_TM_CUT),
    ),
    ## These two MC things depend on the specific choice of probes
    probeMatches  = cms.InputTag("muMcMatch"),
    #allProbes     = cms.InputTag("calProbes"),   # NO 'unbias' efficiency on skims
)

process.allTPHistos = cms.Sequence(
    process.histoMuFromTk +
    process.histoMuFromCal 
)

##    ____       _   _     
##   |  _ \ __ _| |_| |__  
##   | |_) / _` | __| '_ \ 
##   |  __/ (_| | |_| | | |
##   |_|   \__,_|\__|_| |_|
##                         
process.tagAndProbe = cms.Path( 
    process.allProbes *
    process.tagMuons *
    process.tkPassingProbes *
    process.allTPPairs   *
    process.allMcMatches * 
    process.allTPHistos
)

##     ___        _               _     _   _ _     _            
##    / _ \ _   _| |_ _ __  _   _| |_  | | | (_)___| |_ ___  ___ 
##   | | | | | | | __| '_ \| | | | __| | |_| | / __| __/ _ \/ __|
##   | |_| | |_| | |_| |_) | |_| | |_  |  _  | \__ \ || (_) \__ \
##    \___/ \__,_|\__| .__/ \__,_|\__| |_| |_|_|___/\__\___/|___/
##                   |_|                                         
##   
process.TFileService = cms.Service("TFileService", fileName = cms.string("tnpJPsi_MuonID.root"))

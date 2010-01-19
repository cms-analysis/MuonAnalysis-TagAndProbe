import FWCore.ParameterSet.Config as cms

process = cms.Process("TagProbe")

process.load('FWCore.MessageService.MessageLogger_cfi')
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.MessageLogger.cerr.FwkReport.reportEvery = 100

massRange    = ( 2.8, 3.5 )
massRangeSta = ( 2.0, 5.0 )

process.source = cms.Source("PoolSource", 
    fileNames = cms.untracked.vstring(
        #'file:/afs/cern.ch/user/g/gpetrucc/scratch0/huntForRedOctober/CMSSW_3_1_3/src/JPsiMuMu_Skim.root'
        'file:/afs/cern.ch/user/g/gpetrucc/scratch0/tnp/skimJPsi-34X.root'
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
process.GlobalTag.globaltag = cms.string('STARTUP3X_V15::All')

##    ____                   _____               _      ____            _               
##   | __ )  __ _ _ __ ___  |_   _| __ __ _  ___| | __ |  _ \ _ __ ___ | |__   ___  ___ 
##   |  _ \ / _` | '__/ _ \   | || '__/ _` |/ __| |/ / | |_) | '__/ _ \| '_ \ / _ \/ __|
##   | |_) | (_| | | |  __/   | || | | (_| | (__|   <  |  __/| | | (_) | |_) |  __/\__ \
##   |____/ \__,_|_|  \___|   |_||_|  \__,_|\___|_|\_\ |_|   |_|  \___/|_.__/ \___||___/
##                                                                                      
##   
TRACK_CUTS="track.numberOfValidHits >= 10 && track.normalizedChi2 < 5 && abs(track.d0) < 2 && abs(track.dz) < 30"
PT_ETA_CUTS="((pt > 3 && abs(eta) < 1.3) || (pt > 1.5 && 1.3 <= abs(eta) < 2.4))"
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

##    __  __                     ____            _                                 _   _____               
##   |  \/  |_   _  ___  _ __   |  _ \ _ __ ___ | |__   ___  ___    __ _ _ __   __| | |_   _|_ _  __ _ ___ 
##   | |\/| | | | |/ _ \| '_ \  | |_) | '__/ _ \| '_ \ / _ \/ __|  / _` | '_ \ / _` |   | |/ _` |/ _` / __|
##   | |  | | |_| | (_) | | | | |  __/| | | (_) | |_) |  __/\__ \ | (_| | | | | (_| |   | | (_| | (_| \__ \
##   |_|  |_|\__,_|\___/|_| |_| |_|   |_|  \___/|_.__/ \___||___/  \__,_|_| |_|\__,_|   |_|\__,_|\__, |___/
##                                                                                               |___/     
##   
PASS_HLT = "!triggerObjectMatchesByPath('%s').empty()" % ("HLT_Mu3",);
process.tagMuons = cms.EDFilter("PATMuonRefSelector",
    src = cms.InputTag("patMuons"),
    cut = cms.string("isGlobalMuon && pt > 3 && abs(eta) < 2.4 && " + TRACK_CUTS + " && " +PASS_HLT ), 
    filter = cms.bool(True),
)

process.calProbes = cms.EDFilter("PATMuonRefSelector",
    src = cms.InputTag("patMuons"),
#   cut = cms.string("isCaloMuon && pt > 3 && abs(eta) < 2.4 && " + TRACK_CUTS), 
    cut = cms.string("isCaloMuon && "+PT_ETA_CUTS+" && "+TRACK_CUTS),
)
process.glbMuons = cms.EDFilter("PATMuonRefSelector",
    src = cms.InputTag("patMuons"),
    cut = cms.string("isGlobalMuon && "+TRACK_CUTS), 
)
process.glbProbes = cms.EDFilter("PATMuonRefSelector",
    src = cms.InputTag("patMuons"), # can't use glbMuons as source, as RefSelectors can't be chained :-/
    cut = cms.string("isGlobalMuon && pt > 3 && abs(eta) < 2.1 && " + TRACK_CUTS), # 2.1, as we want to use it for trigger!
)
process.trkExMuons = cms.EDFilter("PATMuonRefSelector",
    src = cms.InputTag("patMuons"),
    cut = cms.string("!isGlobalMuon && isTrackerMuon && (isGood('TMLastStationOptimizedLowPtLoose') || isGood('TM2DCompatibilityTight')) && "+TRACK_CUTS),
)


process.allTagsAndProbes = cms.Sequence(
    process.tagMuons +
    process.betterTracks * process.tkTracks * process.tkProbes +
    process.staTracks * process.staProbes +
    process.calProbes +
    process.glbMuons * process.glbProbes +
    process.trkExMuons 
)

##    ____               _               ____            _                   __  __         ___ ____  
##   |  _ \ __ _ ___ ___(_)_ __   __ _  |  _ \ _ __ ___ | |__   ___  ___ _  |  \/  |_   _  |_ _|  _ \ 
##   | |_) / _` / __/ __| | '_ \ / _` | | |_) | '__/ _ \| '_ \ / _ \/ __(_) | |\/| | | | |  | || | | |
##   |  __/ (_| \__ \__ \ | | | | (_| | |  __/| | | (_) | |_) |  __/\__ \_  | |  | | |_| |  | || |_| |
##   |_|   \__,_|___/___/_|_| |_|\__, | |_|   |_|  \___/|_.__/ \___||___(_) |_|  |_|\__,_| |___|____/ 
##                               |___/                                                                
##   
process.tkToGlbMatch = cms.EDProducer("MatcherUsingTracks",
    src     = cms.InputTag("tkTracks"), # all tracks are available for matching
    matched = cms.InputTag("glbMuons"), # to all global muons
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
process.muToGlbMatch = process.tkToGlbMatch.clone(
    src     = cms.InputTag("patMuons"), # again, start with all  muons
    matched = cms.InputTag("glbMuons"), # and match to global
)
process.tkPassingGlb = cms.EDProducer("MatchedCandidateSelector",
    src   = cms.InputTag("tkProbes"),
    match = cms.InputTag("tkToGlbMatch"),
)
process.calPassingGlb = cms.EDProducer("MatchedCandidateSelector",
    src   = cms.InputTag("calProbes"),
    match = cms.InputTag("muToGlbMatch"),
)
process.muToTrkExMatch  = process.muToGlbMatch.clone(matched = 'trkExMuons')
process.tkToTrkExMatch  = process.tkToGlbMatch.clone(matched = 'trkExMuons')
process.tkPassingTrkEx  = process.tkPassingGlb.clone( match = 'tkToTrkExMatch')
process.calPassingTrkEx = process.calPassingGlb.clone(match = 'muToTrkExMatch')
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
    maxDeltaEta      = cms.double(0.2),  # small in eta, which is more precise
    maxDeltaLocalPos = cms.double(100),
    maxDeltaPtRel    = cms.double(3),
    sortBy           = cms.string("deltaR"),
)
process.staPassingTk = cms.EDProducer("MatchedCandidateSelector",
    src   = cms.InputTag("staProbes"),
    match = cms.InputTag("staToTkMatch"),
)

process.allPassingProbes = cms.Sequence(
    process.tkToGlbMatch * process.tkPassingGlb +
    process.muToGlbMatch * process.calPassingGlb +
    process.tkToTrkExMatch * process.tkPassingTrkEx +
    process.muToTrkExMatch * process.calPassingTrkEx +
    process.staToTkMatch   * process.staPassingTk 
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
    cut   = cms.string("%f < mass < %f" % massRange),
)
process.tpGlbCal = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("tagMuons@+ calProbes@-"), # charge coniugate states are implied
    cut   = cms.string("%f < mass < %f" % massRange),
)
process.tpGlbSta = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("tagMuons@+ staProbes@-"), # charge coniugate states are implied
    cut   = cms.string("%f < mass < %f" % massRangeSta),
)
process.tpGlbGlb = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("tagMuons@+ glbProbes@-"), # charge coniugate states are implied
    cut   = cms.string("%f < mass < %f" % massRange),
)
process.allTPPairs = cms.Sequence(process.tpGlbTk + process.tpGlbCal + process.tpGlbSta + process.tpGlbGlb)


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
process.staMcMatch = process.muMcMatch.clone(src = "staTracks", distMin = 0.6)

process.allMcMatches = cms.Sequence(process.muMcMatch + process.tkMcMatch + process.staMcMatch)


##    _____           _       _ ____            _            _   _ _____            _      
##   |_   _|_ _  __ _( )_ __ ( )  _ \ _ __ ___ | |__   ___  | \ | |_   _|   _ _ __ | | ___ 
##     | |/ _` |/ _` |/| '_ \|/| |_) | '__/ _ \| '_ \ / _ \ |  \| | | || | | | '_ \| |/ _ \
##     | | (_| | (_| | | | | | |  __/| | | (_) | |_) |  __/ | |\  | | || |_| | |_) | |  __/
##     |_|\__,_|\__, | |_| |_| |_|   |_|  \___/|_.__/ \___| |_| \_| |_| \__,_| .__/|_|\___|
##              |___/                                                        |_|           
##
recoCommonStuff = cms.PSet(
    variables = cms.PSet(
        eta = cms.string("eta"),
        pt  = cms.string("pt"),
        phi  = cms.string("phi"),
    )
)
mcTruthCommonStuff = cms.PSet(
    isMC = cms.bool(True),
    makeMCUnbiasTree = cms.bool(True),
    checkMotherInUnbiasEff = cms.bool(True),
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
        passingGlb   = cms.InputTag("tkPassingGlb"),
        passingTrkEx = cms.InputTag("tkPassingTrkEx"),
    ),
    ## These two MC things depend on the specific choice of probes
    probeMatches  = cms.InputTag("tkMcMatch"),
    allProbes     = cms.InputTag("tkProbes"),
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
        passingGlb   = cms.InputTag("calPassingGlb"),
        passingTrkEx = cms.InputTag("calPassingTrkEx"),
    ),
    ## These two MC things depend on the specific choice of probes
    probeMatches  = cms.InputTag("muMcMatch"),
    allProbes     = cms.InputTag("calProbes"),
)

#####
## Tk from Sta
process.histoTkFromSta = cms.EDAnalyzer("TagProbeFitTreeProducer",
    ## pick the defaults
    recoCommonStuff, mcTruthCommonStuff,
    # choice of tag and probe pairs, and arbitration
    tagProbePairs = cms.InputTag("tpGlbSta"),
    arbitration   = cms.string("OneProbe"),
    # choice of what defines a 'passing' probe
    flags = cms.PSet(passing = cms.InputTag("staPassingTk")),
    ## These two MC things depend on the specific choice of probes
    probeMatches  = cms.InputTag("staMcMatch"),
    allProbes     = cms.InputTag("staProbes"),
)

#####
## HLT from Glb
process.histoHltFromGlb = cms.EDAnalyzer("TagProbeFitTreeProducer",
    ## pick the defaults
    recoCommonStuff, mcTruthCommonStuff,
    # choice of tag and probe pairs, and arbitration
    tagProbePairs = cms.InputTag("tpGlbGlb"),
    arbitration   = cms.string("OneProbe"),
    # choice of what defines a 'passing' probe
    flags = cms.PSet(
        L1MuOpen = cms.string("!triggerObjectMatchesByPath('HLT_L1MuOpen').empty()"),
        Mu3 = cms.string("!triggerObjectMatchesByPath('HLT_Mu3').empty()"),
        Mu5 = cms.string("!triggerObjectMatchesByPath('HLT_Mu5').empty()"),
        DoubleMu0 = cms.string("!triggerObjectMatchesByPath('HLT_DoubleMu0').empty()"),
        DoubleMu3 = cms.string("!triggerObjectMatchesByPath('HLT_DoubleMu3').empty()"),
    ),
    ## These two MC things depend on the specific choice of probes
    probeMatches  = cms.InputTag("muMcMatch"),
    allProbes     = cms.InputTag("glbProbes"),
)


process.allTPHistos = cms.Sequence(
    process.histoMuFromTk +
    process.histoMuFromCal +
    process.histoTkFromSta +
    process.histoHltFromGlb 
)

##    ____       _   _     
##   |  _ \ __ _| |_| |__  
##   | |_) / _` | __| '_ \ 
##   |  __/ (_| | |_| | | |
##   |_|   \__,_|\__|_| |_|
##                         
process.tagAndProbe = cms.Path( 
    process.allTagsAndProbes *
    process.allPassingProbes *
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
process.TFileService = cms.Service("TFileService", fileName = cms.string("tnpJPsi.root"))


##    ____             _      ___         _____ _                
##   | __ )  __ _  ___| | __ |_ _|_ __   |_   _(_)_ __ ___   ___ 
##   |  _ \ / _` |/ __| |/ /  | || '_ \    | | | | '_ ` _ \ / _ \
##   | |_) | (_| | (__|   <   | || | | |   | | | | | | | | |  __/
##   |____/ \__,_|\___|_|\_\ |___|_| |_|   |_| |_|_| |_| |_|\___|
##                                                               
#### To run on 3.1.X Skims
#### change output file and test input file
# process.source.fileNames = [ 'file:/afs/cern.ch/user/g/gpetrucc/scratch0/huntForRedOctober/CMSSW_3_1_3/src/JPsiMuMu_Skim.root' ]
# process.TFileService.fileName = 'tnpJPsi31X.root'
#### we had only HLT_Mu3
# process.histoHltFromGlb.flags = cms.PSet(Mu3 = cms.string("!triggerObjectMatchesByPath('HLT_Mu3').empty()"))

#
##    ____       _                 
##   |  _ \  ___| |__  _   _  __ _ 
##   | | | |/ _ \ '_ \| | | |/ _` |
##   | |_| |  __/ |_) | |_| | (_| |
##   |____/ \___|_.__/ \__,_|\__, |
##                           |___/ 
##

### limit events
#process.maxEvents.input = 200



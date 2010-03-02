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
PASSING_TRK_CUT = "!("+PASSING_GLB_CUT+") && isTrackerMuon && muonID('TMLastStationAngTight')";

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

process.tagMuons = cms.EDFilter("PATMuonRefSelector",
    src = cms.InputTag("patMuons"),
    cut = cms.string(TAG_CUTS), 
    filter = cms.bool(True),
)

process.probeMuons = cms.EDFilter("PATMuonRefSelector",
    src = cms.InputTag("patMuons"),
    cut = cms.string("track.isNonnull && "+TRACK_CUTS), # pick everything here, select Glb/Trk later
)

process.tp = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("tagMuons@+ probeMuons@-"), # charge coniugate states are implied
    cut   = cms.string("%f < mass < %f" % MASS_RANGE),
)

process.muMcMatch = cms.EDFilter("MCTruthDeltaRMatcherNew",
    pdgId = cms.vint32(13),
    src = cms.InputTag("patMuons"),
    distMin = cms.double(0.3),
    matched = cms.InputTag("genMuons")
)

process.histo = cms.EDAnalyzer("TagProbeFitTreeProducer",
    variables = cms.PSet(
        pt  = cms.string("pt"),
        p   = cms.string("p"),
        eta = cms.string("eta"),
        phi = cms.string("phi"),
    ),
    # choice of tag and probe pairs, and arbitration
    tagProbePairs = cms.InputTag("tp"),
    arbitration   = cms.string("OneProbe"),
    # choice of what defines a 'passing' probe
    flags = cms.PSet(
        Glb = cms.string(PASSING_GLB_CUT),
        Trk = cms.string(PASSING_TRK_CUT),
        HLTMu3     = cms.string("!triggerObjectMatchesByFilter('hltSingleMu3L3Filtered3').empty()"),
        L1DiMuOpen = cms.string("!triggerObjectMatchesByFilter('hltDoubleMuLevel1PathL1OpenFiltered').empty()"),
    ),
    isMC = cms.bool(True),
    tagMatches = cms.InputTag("muMcMatch"),
    probeMatches  = cms.InputTag("muMcMatch"),
    motherPdgId = cms.int32(443),
    makeMCUnbiasTree = cms.bool(True),
    checkMotherInUnbiasEff = cms.bool(True),
    allProbes = cms.InputTag("probeMuons"),
)
           
process.tagAndProbe = cms.Path( 
    process.tagMuons *
    process.probeMuons *
    process.tp *
    process.muMcMatch   *
    process.histo 
)

process.TFileService = cms.Service("TFileService", fileName = cms.string("tnpJPsi_Trigger.root"))

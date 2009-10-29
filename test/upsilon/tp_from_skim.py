import FWCore.ParameterSet.Config as cms

process = cms.Process("TagProbe")

process.load('FWCore.MessageService.MessageLogger_cfi')
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.MessageLogger.cerr.FwkReport.reportEvery = 100

massRange    = ( 8.4, 11.4 ); 
massRangeSta = (5, 14)

process.source = cms.Source("PoolSource", 
    fileNames = cms.untracked.vstring(
    'file:/tmp/jjhollar/skimUpsilon_1.root'
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
process.GlobalTag.globaltag = cms.string('MC_31X_V3::All')

##    ____                   _____               _      ____            _               
##   | __ )  __ _ _ __ ___  |_   _| __ __ _  ___| | __ |  _ \ _ __ ___ | |__   ___  ___ 
##   |  _ \ / _` | '__/ _ \   | || '__/ _` |/ __| |/ / | |_) | '__/ _ \| '_ \ / _ \/ __|
##   | |_) | (_| | | |  __/   | || | | (_| | (__|   <  |  __/| | | (_) | |_) |  __/\__ \
##   |____/ \__,_|_|  \___|   |_||_|  \__,_|\___|_|\_\ |_|   |_|  \___/|_.__/ \___||___/
##                                                                                      
##   
process.tkTracks  = cms.EDProducer("ConcreteChargedCandidateProducer", 
    src  = cms.InputTag("goodTracks"),      
    particleType = cms.string("mu+"),
) 
process.tkProbes = cms.EDProducer("CandViewRefSelector",
    src = cms.InputTag("tkTracks"),
    cut = cms.string("pt > 3 && abs(eta) < 2.4"),
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
    cut = cms.string("pt > 3 && abs(eta) < 2.4"),
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
    cut = cms.string("isGlobalMuon && pt > 3 && abs(eta) < 2.4 && " + PASS_HLT ), 
    filter = cms.bool(True),
)

process.calProbes = cms.EDFilter("PATMuonRefSelector",
    src = cms.InputTag("patMuons"),
    cut = cms.string("isCaloMuon && pt > 3 && abs(eta) < 2.4"), 
)
process.glbMuons = cms.EDFilter("PATMuonRefSelector",
    src = cms.InputTag("patMuons"),
    cut = cms.string("isGlobalMuon"), 
)
process.glbProbes = cms.EDFilter("PATMuonRefSelector",
    src = cms.InputTag("patMuons"), # can't use glbMuons as source, as RefSelectors can't be chained :-/
    cut = cms.string("isGlobalMuon && pt > 3 && abs(eta) < 2.1"), # 2.1, as we want to use it for trigger!
)

process.allTagsAndProbes = cms.Sequence(
    process.tagMuons +
    process.tkTracks * process.tkProbes +
    process.staTracks * process.staProbes +
    process.calProbes +
    process.glbMuons * process.glbProbes 
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

##    ____               _               ____            _                   _____     _                       
##   |  _ \ __ _ ___ ___(_)_ __   __ _  |  _ \ _ __ ___ | |__   ___  ___ _  |_   _| __(_) __ _  __ _  ___ _ __ 
##   | |_) / _` / __/ __| | '_ \ / _` | | |_) | '__/ _ \| '_ \ / _ \/ __(_)   | || '__| |/ _` |/ _` |/ _ \ '__|
##   |  __/ (_| \__ \__ \ | | | | (_| | |  __/| | | (_) | |_) |  __/\__ \_    | || |  | | (_| | (_| |  __/ |   
##   |_|   \__,_|___/___/_|_| |_|\__, | |_|   |_|  \___/|_.__/ \___||___(_)   |_||_|  |_|\__, |\__, |\___|_|   
##                               |___/                                                   |___/ |___/           
##   
process.glbPassingHLT = cms.EDFilter("PATMuonRefSelector",
    src = cms.InputTag("patMuons"),
    cut = cms.string(process.glbProbes.cut.value() + " && " + PASS_HLT),
)

##    ____               _               ____            _                   ____        ____   _    ____ 
##   |  _ \ __ _ ___ ___(_)_ __   __ _  |  _ \ _ __ ___ | |__   ___  ___ _  | __ )      |  _ \ / \  / ___|
##   | |_) / _` / __/ __| | '_ \ / _` | | |_) | '__/ _ \| '_ \ / _ \/ __(_) |  _ \ _____| |_) / _ \| |  _ 
##   |  __/ (_| \__ \__ \ | | | | (_| | |  __/| | | (_) | |_) |  __/\__ \_  | |_) |_____|  __/ ___ \ |_| |
##   |_|   \__,_|___/___/_|_| |_|\__, | |_|   |_|  \___/|_.__/ \___||___(_) |____/      |_| /_/   \_\____|
##                               |___/                                                                    
##   

## 1) From the J/Psi + Mu analysis
process.muonsJPsiPlusMu = cms.EDProducer("MuonSelectorJPsiPlusMu",
    src = cms.InputTag("patMuons")
)
process.tkToMuJPsiPlusMuMatch = process.tkToGlbMatch.clone(matched = "muonsJPsiPlusMu")
process.tkPassingMuJPsiPlusMu = process.tkPassingGlb.clone(match   = "tkToMuJPsiPlusMuMatch")
## 2a) From the J/Psi analysis (Global)
process.muonsJPsiGlb = cms.EDProducer("MuonSelectorJPsi",
    src = cms.InputTag("patMuons"),
    selectGlobalMuons = cms.bool(True)
)
process.tkToMuJPsiGlbMatch = process.tkToGlbMatch.clone(matched = "muonsJPsiGlb")
process.tkPassingMuJPsiGlb = process.tkPassingGlb.clone(match   = "tkToMuJPsiGlbMatch")

## 2b) From the J/Psi analysis (Tracker not global)
process.muonsJPsiTkM = cms.EDProducer("MuonSelectorJPsi",
    src = cms.InputTag("patMuons"),
    selectGlobalMuons = cms.bool(True)
)
process.tkToMuJPsiTkMMatch = process.tkToGlbMatch.clone(matched = "muonsJPsiTkM")
process.tkPassingMuJPsiTkM = process.tkPassingGlb.clone(match   = "tkToMuJPsiTkMMatch")

## 2c) From the J/Psi analysis (Global with custom ID cut -> trigger)
process.muJPsiGlbProbes = cms.EDProducer("MuonSelectorJPsi",
    src = cms.InputTag("glbProbes"),
    selectGlobalMuons = cms.bool(True)
)
process.muJPsiGlbPassingHLT = cms.EDProducer("MuonSelectorJPsi",
    src = cms.InputTag("glbPassingHLT"),
    selectGlobalMuons = cms.bool(True)
)


## 3) From the Exclusive B analysis: a very simple cut, just the OR of global and tracker muon
process.muonsBExcl = cms.EDFilter("PATMuonRefSelector",
    src = cms.InputTag("patMuons"),
    cut = cms.string("isGlobalMuon || isTrackerMuon"),
)
process.tkToMuBExclMatch = process.tkToGlbMatch.clone(matched = "muonsBExcl")
process.tkPassingMuBExcl = process.tkPassingGlb.clone(match   = "tkToMuBExclMatch")

process.PAGMuonIDs = cms.Sequence(
    process.muonsJPsiPlusMu * process.tkToMuJPsiPlusMuMatch * process.tkPassingMuJPsiPlusMu +
    process.muonsJPsiGlb * process.tkToMuJPsiGlbMatch * process.tkPassingMuJPsiGlb +
    process.muJPsiGlbProbes * process.muJPsiGlbPassingHLT +
    process.muonsJPsiTkM * process.tkToMuJPsiTkMMatch * process.tkPassingMuJPsiTkM +
    process.muonsBExcl * process.tkToMuBExclMatch * process.tkPassingMuBExcl 
)

process.allPassingProbes = cms.Sequence(
    process.tkToGlbMatch * process.tkPassingGlb +
    process.muToGlbMatch * process.calPassingGlb +
    process.staToTkMatch * process.staPassingTk +
    process.glbPassingHLT +
    process.PAGMuonIDs 
)



##    __  __       _          _____                 ____            _                                      
##   |  \/  | __ _| | _____  |_   _|_ _  __ _   _  |  _ \ _ __ ___ | |__   ___   _ __ ___   __ _ _ __  ___ 
##   | |\/| |/ _` | |/ / _ \   | |/ _` |/ _` |_| |_| |_) | '__/ _ \| '_ \ / _ \ | '_ ` _ \ / _` | '_ \/ __|
##   | |  | | (_| |   <  __/   | | (_| | (_| |_   _|  __/| | | (_) | |_) |  __/ | | | | | | (_| | |_) \__ \
##   |_|  |_|\__,_|_|\_\___|   |_|\__,_|\__, | |_| |_|   |_|  \___/|_.__/ \___| |_| |_| |_|\__,_| .__/|___/
##                                      |___/                                                   |_|        
tagProbeTemplate = cms.EDProducer("TagProbeProducer",
    TagCollection = cms.InputTag("tagMuons"),
    MassMinCut = cms.untracked.double(massRange[0]),
    MassMaxCut = cms.untracked.double(massRange[1]),
)

process.tagProbeGlbFromTk  = tagProbeTemplate.clone(  ProbeCollection = cms.InputTag("tkProbes")  )
process.tagProbeGlbFromCal = tagProbeTemplate.clone(  ProbeCollection = cms.InputTag("calProbes") )
process.tagProbeTkFromSta  = tagProbeTemplate.clone(  ProbeCollection = cms.InputTag("staProbes"),
                                                      MassMinCut = cms.untracked.double(massRangeSta[0]),
                                                      MassMaxCut = cms.untracked.double(massRangeSta[1]), )
process.tagProbeHltFromGlb = tagProbeTemplate.clone(  ProbeCollection = cms.InputTag("glbProbes") )

### Note that we do NOT need extra TagProbeProducers, as the PAG-specific muon IDs just define new *passing* probes.
### Except for the Global->HLT
process.tagProbeHltFromJPsiGlb = tagProbeTemplate.clone(  ProbeCollection = cms.InputTag("muJPsiGlbProbes") )


process.allTagProbeMaps = cms.Sequence(
    process.tagProbeGlbFromTk +
    process.tagProbeGlbFromCal +
    process.tagProbeTkFromSta +
    process.tagProbeHltFromGlb +
    process.tagProbeHltFromJPsiGlb
)

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
process.TPEdm = cms.EDProducer("TagProbeEDMNtuple",
    tagProbeType = cms.untracked.string( "Muon" ),
    genParticlesTag = cms.untracked.InputTag("genMuons"),
    mcParticles = cms.untracked.vint32(13),
    mcParents   = cms.untracked.vint32( 0),
    
    ## Tag & Probe Muon Candidate Collections
    ## all this will be filled once for each T&P measurement
    tagCandTags = cms.untracked.VInputTag(),
    passProbeCandTags = cms.untracked.VInputTag(),
    allProbeCandTags  = cms.untracked.VInputTag(),
    tagTruthMatchMapTags       = cms.untracked.VInputTag(),
    passProbeTruthMatchMapTags = cms.untracked.VInputTag(),
    allProbeTruthMatchMapTags  = cms.untracked.VInputTag(),
    tagProbeMapTags = cms.untracked.VInputTag(),
    BestProbeCriteria = cms.untracked.vstring(),
    BestProbeInvMass  = cms.untracked.vdouble(),

    ## Ignore the following, we don't care
    triggerEventTag = cms.untracked.InputTag( "hltTriggerSummaryAOD" ),
    hltL1Tag        = cms.untracked.InputTag( "hltSingleMuIsoLevel1Seed" ),
    hltTag          = cms.untracked.InputTag( "hltSingleMuIsoL3IsoFiltered" ),
    triggerDelRMatch     = cms.untracked.double( 0.15 ),
    triggerDelPtRelMatch = cms.untracked.double( 0.15 )
)

##    _____           _       _ ____            _            _   _ _     _        
##   |_   _|_ _  __ _( )_ __ ( )  _ \ _ __ ___ | |__   ___  | | | (_)___| |_ ___  
##     | |/ _` |/ _` |/| '_ \|/| |_) | '__/ _ \| '_ \ / _ \ | |_| | / __| __/ _ \ 
##     | | (_| | (_| | | | | | |  __/| | | (_) | |_) |  __/ |  _  | \__ \ || (_) |
##     |_|\__,_|\__, | |_| |_| |_|   |_|  \___/|_.__/ \___| |_| |_|_|___/\__\___/ 
##              |___/                                                             
MakeHisto = cms.EDAnalyzer("TagProbeEDMAnalysis",      
      # AFAIK these are the only meaningful parameters for now that we just make the Histo
      # --------------------------------------------
      Mode = cms.untracked.string("Write"),
      FitFileName = cms.untracked.string( "fit_output.root" ),

      TagProbeType = cms.untracked.int32(0),

      MCTruthParentId = cms.untracked.int32(443),
      Weight = cms.untracked.double(1.0),

      CalculateEffSideBand = cms.untracked.bool( True ), ## Calculate and store effs using SB
      CalculateEffFitter   = cms.untracked.bool( True ), ## Calculate and store effs from Roofit
      CalculateEffTruth    = cms.untracked.bool( True ), ## Calculate and store true effs
      UnbinnedFit          = cms.untracked.bool( True ),
      Do2DFit              = cms.untracked.bool( True ),

      NumBinsMass         = cms.untracked.int32( 20 ),
      MassLow             = cms.untracked.double( massRange[0] ),
      MassHigh            = cms.untracked.double( massRange[1] ),

      NameVar1             = cms.untracked.string( "pt" ),
      Var1BinBoundaries   = cms.untracked.vdouble( 3, 4.5, 6, 10, 20 ),
      NameVar2             = cms.untracked.string( "eta" ),
      Var2BinBoundaries   = cms.untracked.vdouble( -2.4,-1.2,-0.7,0.0,0.7,1.2,2.4),

      # All the following is useless now
      # --------------------------------------------
      GaussLineShape = cms.untracked.PSet(
        GaussMean        = cms.untracked.vdouble( 3.09, 2.9,  3.1 ),
        GaussSigma       = cms.untracked.vdouble( 0.03, 0.01, 0.05 )
      ),

      CMSBkgLineShape = cms.untracked.PSet(
        CMSBkgAlpha           = cms.untracked.vdouble( 124, 0, 1000 ),
        CMSBkgBeta            = cms.untracked.vdouble( -0.028,-1.0,1.0 ),
        CMSBkgPeak            = cms.untracked.vdouble( 91.1876 ),
        CMSBkgGamma           = cms.untracked.vdouble( 0.0379,0.0,0.5 )
      ),

      Efficiency        = cms.untracked.vdouble( 0.99,0.5,1.0 ),    
      NumSignal         = cms.untracked.vdouble( 4000.0,-10.0,30000.0 ),    
      NumBkgPass        = cms.untracked.vdouble( 4000.0,0.0,10000.0 ),    
      NumBkgFail        = cms.untracked.vdouble( 1000.0,-10.0,7000.0 ),    

      SBSPeak            = cms.untracked.double( 3.1 ),   ## Mass peak
      SBSStanDev         = cms.untracked.double( 2 )        ## SD from peak for subtraction
)

##    __  __       _                          _       _____      ____                        _       _      
##   |  \/  | __ _| | _____    ___  __ _  ___| |__   |_   _| __ |  _ \   _ __ ___   ___   __| |_   _| | ___ 
##   | |\/| |/ _` | |/ / _ \  / _ \/ _` |/ __| '_ \    | || '_ \| |_) | | '_ ` _ \ / _ \ / _` | | | | |/ _ \
##   | |  | | (_| |   <  __/ |  __/ (_| | (__| | | |   | || | | |  __/  | | | | | | (_) | (_| | |_| | |  __/
##   |_|  |_|\__,_|_|\_\___|  \___|\__,_|\___|_| |_|   |_||_| |_|_|     |_| |_| |_|\___/ \__,_|\__,_|_|\___|
##                                                                                                          

#####
## Mu from Tk
process.TPEdm.tagCandTags       += [ cms.InputTag("tagMuons") ]
process.TPEdm.allProbeCandTags  += [ cms.InputTag("tkProbes")   ]
process.TPEdm.passProbeCandTags += [ cms.InputTag("tkPassingGlb") ]
process.TPEdm.tagProbeMapTags   += [ cms.InputTag("tagProbeGlbFromTk") ]
process.TPEdm.tagTruthMatchMapTags       += [ cms.InputTag("muMcMatch") ]
process.TPEdm.passProbeTruthMatchMapTags += [ cms.InputTag("tkMcMatch") ]
process.TPEdm.allProbeTruthMatchMapTags  += [ cms.InputTag("tkMcMatch") ]
process.TPEdm.BestProbeCriteria += [ "OneProbe" ]
process.TPEdm.BestProbeInvMass  += [ 3.1 ]
process.fitGlbFromTk = MakeHisto.clone( 
    TagProbeType = cms.untracked.int32(0),
    FitFileName = cms.untracked.string( "histo_output_GlbFromTk.root"),
)

#####
## Mu from Cal
process.TPEdm.tagCandTags       += [ cms.InputTag("tagMuons") ]
process.TPEdm.allProbeCandTags  += [ cms.InputTag("calProbes")   ]
process.TPEdm.passProbeCandTags += [ cms.InputTag("calPassingGlb") ]
process.TPEdm.tagProbeMapTags   += [ cms.InputTag("tagProbeGlbFromCal") ]
process.TPEdm.tagTruthMatchMapTags       += [ cms.InputTag("muMcMatch") ]
process.TPEdm.passProbeTruthMatchMapTags += [ cms.InputTag("muMcMatch") ]
process.TPEdm.allProbeTruthMatchMapTags  += [ cms.InputTag("muMcMatch") ]
process.TPEdm.BestProbeCriteria += [ "OneProbe" ]
process.TPEdm.BestProbeInvMass  += [ 3.1 ]
process.fitGlbFromCal = MakeHisto.clone( 
    TagProbeType = cms.untracked.int32(1),
    FitFileName = cms.untracked.string( "histo_output_GlbFromCal.root"),
)

#####
## Tk from Sta
process.TPEdm.tagCandTags       += [ cms.InputTag("tagMuons") ]
process.TPEdm.allProbeCandTags  += [ cms.InputTag("staProbes")   ]
process.TPEdm.passProbeCandTags += [ cms.InputTag("staPassingTk") ]
process.TPEdm.tagProbeMapTags   += [ cms.InputTag("tagProbeTkFromSta") ]
process.TPEdm.tagTruthMatchMapTags       += [ cms.InputTag("muMcMatch") ]
process.TPEdm.passProbeTruthMatchMapTags += [ cms.InputTag("staMcMatch") ]
process.TPEdm.allProbeTruthMatchMapTags  += [ cms.InputTag("staMcMatch") ]
process.TPEdm.BestProbeCriteria += [ "OneProbe" ]
process.TPEdm.BestProbeInvMass  += [ 3.1 ]
process.fitTkFromSta = MakeHisto.clone( 
    TagProbeType = cms.untracked.int32(2),
    FitFileName = cms.untracked.string( "histo_output_TkFromSta.root"),
    ## need to override mass range
    MassLow  = massRangeSta[0],   
    MassHigh = massRangeSta[1],
)

#####
## HLT from Glb
process.TPEdm.tagCandTags       += [ cms.InputTag("tagMuons") ]
process.TPEdm.allProbeCandTags  += [ cms.InputTag("glbProbes")   ]
process.TPEdm.passProbeCandTags += [ cms.InputTag("glbPassingHLT") ]
process.TPEdm.tagProbeMapTags   += [ cms.InputTag("tagProbeHltFromGlb") ]
process.TPEdm.tagTruthMatchMapTags       += [ cms.InputTag("muMcMatch") ]
process.TPEdm.passProbeTruthMatchMapTags += [ cms.InputTag("muMcMatch") ]
process.TPEdm.allProbeTruthMatchMapTags  += [ cms.InputTag("muMcMatch") ]
process.TPEdm.BestProbeCriteria += [ "OneProbe" ]
process.TPEdm.BestProbeInvMass  += [ 3.1 ]
process.fitHltFromGlb = MakeHisto.clone( 
    TagProbeType = cms.untracked.int32(3),
    FitFileName = cms.untracked.string( "histo_output_HltFromGlb.root"),
    ## need to override eta bins
    Var2BinBoundaries   = cms.untracked.vdouble( -2.1,-1.2,-0.7,0.0,0.7,1.2,2.1 ),
)

process.allTPHistos = cms.Sequence(
    process.TPEdm +
    process.fitGlbFromTk +
    process.fitGlbFromCal +
    process.fitTkFromSta +
    process.fitHltFromGlb 
)

## PAG-specific muon IDs
process.TPEdm.tagCandTags       += [ cms.InputTag("tagMuons") ]
process.TPEdm.allProbeCandTags  += [ cms.InputTag("tkProbes")   ]
process.TPEdm.passProbeCandTags += [ cms.InputTag("tkPassingMuJPsiGlb") ]
process.TPEdm.tagProbeMapTags   += [ cms.InputTag("tagProbeGlbFromTk") ]
process.TPEdm.tagTruthMatchMapTags       += [ cms.InputTag("muMcMatch") ]
process.TPEdm.passProbeTruthMatchMapTags += [ cms.InputTag("tkMcMatch") ]
process.TPEdm.allProbeTruthMatchMapTags  += [ cms.InputTag("tkMcMatch") ]
process.TPEdm.BestProbeCriteria += [ "OneProbe" ]
process.TPEdm.BestProbeInvMass  += [ 3.1 ]
process.fitMuFromTkJPsiGlb = MakeHisto.clone( 
    TagProbeType = cms.untracked.int32(4),
    FitFileName = cms.untracked.string( "histo_output_MuFromTkJPsiGlb.root"),
)

process.TPEdm.tagCandTags       += [ cms.InputTag("tagMuons") ]
process.TPEdm.allProbeCandTags  += [ cms.InputTag("tkProbes")   ]
process.TPEdm.passProbeCandTags += [ cms.InputTag("tkPassingMuJPsiTkM") ]
process.TPEdm.tagProbeMapTags   += [ cms.InputTag("tagProbeGlbFromTk") ]
process.TPEdm.tagTruthMatchMapTags       += [ cms.InputTag("muMcMatch") ]
process.TPEdm.passProbeTruthMatchMapTags += [ cms.InputTag("tkMcMatch") ]
process.TPEdm.allProbeTruthMatchMapTags  += [ cms.InputTag("tkMcMatch") ]
process.TPEdm.BestProbeCriteria += [ "OneProbe" ]
process.TPEdm.BestProbeInvMass  += [ 3.1 ]
process.fitMuFromTkJPsiTkM = MakeHisto.clone( 
    TagProbeType = cms.untracked.int32(4),
    FitFileName = cms.untracked.string( "histo_output_MuFromTkJPsiTkM.root"),
)

process.TPEdm.tagCandTags       += [ cms.InputTag("tagMuons") ]
process.TPEdm.allProbeCandTags  += [ cms.InputTag("muJPsiGlbProbes")   ]
process.TPEdm.passProbeCandTags += [ cms.InputTag("muJPsiGlbPassingHLT") ]
process.TPEdm.tagProbeMapTags   += [ cms.InputTag("tagProbeHltFromJPsiGlb") ]
process.TPEdm.tagTruthMatchMapTags       += [ cms.InputTag("muMcMatch") ]
process.TPEdm.passProbeTruthMatchMapTags += [ cms.InputTag("muMcMatch") ]
process.TPEdm.allProbeTruthMatchMapTags  += [ cms.InputTag("muMcMatch") ]
process.TPEdm.BestProbeCriteria += [ "OneProbe" ]
process.TPEdm.BestProbeInvMass  += [ 3.1 ]
process.fitHltFromJPsiGlb = MakeHisto.clone( 
    TagProbeType = cms.untracked.int32(5),
    FitFileName = cms.untracked.string( "histo_output_HltFromJPsiGlb.root"),
    ## need to override eta bins
    Var2BinBoundaries   = cms.untracked.vdouble( -2.1,-1.2,-0.7,0.0,0.7,1.2,2.1 ),
)

process.TPEdm.tagCandTags       += [ cms.InputTag("tagMuons") ]
process.TPEdm.allProbeCandTags  += [ cms.InputTag("tkProbes")   ]
process.TPEdm.passProbeCandTags += [ cms.InputTag("tkPassingMuJPsiPlusMu") ]
process.TPEdm.tagProbeMapTags   += [ cms.InputTag("tagProbeGlbFromTk") ]
process.TPEdm.tagTruthMatchMapTags       += [ cms.InputTag("muMcMatch") ]
process.TPEdm.passProbeTruthMatchMapTags += [ cms.InputTag("tkMcMatch") ]
process.TPEdm.allProbeTruthMatchMapTags  += [ cms.InputTag("tkMcMatch") ]
process.TPEdm.BestProbeCriteria += [ "OneProbe" ]
process.TPEdm.BestProbeInvMass  += [ 3.1 ]
process.fitMuFromTkJPsiPlusMu = MakeHisto.clone( 
    TagProbeType = cms.untracked.int32(6),
    FitFileName = cms.untracked.string( "histo_output_MuFromTkJPsiPlusMu.root"),
)

process.TPEdm.tagCandTags       += [ cms.InputTag("tagMuons") ]
process.TPEdm.allProbeCandTags  += [ cms.InputTag("tkProbes")   ]
process.TPEdm.passProbeCandTags += [ cms.InputTag("tkPassingMuBExcl") ]
process.TPEdm.tagProbeMapTags   += [ cms.InputTag("tagProbeGlbFromTk") ]
process.TPEdm.tagTruthMatchMapTags       += [ cms.InputTag("muMcMatch") ]
process.TPEdm.passProbeTruthMatchMapTags += [ cms.InputTag("tkMcMatch") ]
process.TPEdm.allProbeTruthMatchMapTags  += [ cms.InputTag("tkMcMatch") ]
process.TPEdm.BestProbeCriteria += [ "OneProbe" ]
process.TPEdm.BestProbeInvMass  += [ 3.1 ]
process.fitMuFromTkBExcl = MakeHisto.clone( 
    TagProbeType = cms.untracked.int32(7),
    FitFileName = cms.untracked.string( "histo_output_MuFromTkBExcl.root"),
)

process.allTPHistos += cms.Sequence(
    process.fitMuFromTkJPsiGlb +
    process.fitHltFromJPsiGlb +
    process.fitMuFromTkJPsiTkM +
    process.fitMuFromTkJPsiPlusMu +
    process.fitMuFromTkBExcl 
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
    process.allTagProbeMaps * 
    process.allMcMatches * 
    process.allTPHistos
)

##     ___        _               _   
##    / _ \ _   _| |_ _ __  _   _| |_ 
##   | | | | | | | __| '_ \| | | | __|
##   | |_| | |_| | |_| |_) | |_| | |_ 
##    \___/ \__,_|\__| .__/ \__,_|\__|
##                   |_|              
##
## We define this module, but we don't run it unless needed
process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string("/tmp/gpetrucc/tuple.root"),
    outputCommands = cms.untracked.vstring("drop *", "keep _TPEdm_*_*"),
    dropMetaDataForDroppedData = cms.untracked.bool(True),
)

##    ____       _                 
##   |  _ \  ___| |__  _   _  __ _ 
##   | | | |/ _ \ '_ \| | | |/ _` |
##   | |_| |  __/ |_) | |_| | (_| |
##   |____/ \___|_.__/ \__,_|\__, |
##                           |___/ 
##

### limit events
#process.maxEvents.input = 200




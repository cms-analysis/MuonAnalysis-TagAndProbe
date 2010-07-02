import FWCore.ParameterSet.Config as cms

##   __     __         _       _     _                             _    ____      _       
##   \ \   / /_ _ _ __(_) __ _| |__ | | ___  ___    __ _ _ __   __| |  / ___|   _| |_ ___ 
##    \ \ / / _` | '__| |/ _` | '_ \| |/ _ \/ __|  / _` | '_ \ / _` | | |  | | | | __/ __|
##     \ V / (_| | |  | | (_| | |_) | |  __/\__ \ | (_| | | | | (_| | | |__| |_| | |_\__ \
##      \_/ \__,_|_|  |_|\__,_|_.__/|_|\___||___/  \__,_|_| |_|\__,_|  \____\__,_|\__|___/
##                                                                                        
##   

MASS_RANGE = ( 2.6, 3.6 )
MASS_RANGE_STA = ( 2, 5 )
TRACK_CUTS = ("track.numberOfValidHits > 11 && track.hitPattern.pixelLayersWithMeasurement > 1 "+
              "&& track.normalizedChi2 < 4 "+
              "&& abs(track.d0) < 3 && abs(track.dz) < 30")
PT_ETA_CUTS = "(pt > 3 || (abs(eta)>1 && p > 2.6))" ## the enclosing () are very important, because there's an "||"
JPSI_ACCEPTANCE_CUT = ("(       abs(eta) <= 1.3  && pt > 3.3 || "+
                       "  1.3 < abs(eta) <= 2.2  && p >  2.9 || "+
                       "  2.2 < abs(eta) <= 2.4  && pt > 0.8 )  ")

PASSING_GLB_CUT = ("isGlobalMuon && globalTrack.normalizedChi2 < 20  && "+
                   "globalTrack.hitPattern.numberOfValidMuonHits > 0 && "+
                   "muonID('TrackerMuonArbitrated') &&  muonID('TMLastStationAngTight')")
PASSING_TMI_CUT = "isTrackerMuon && muonID('TMLastStationAngTight')";
PASSING_TM_CUT  = "!("+PASSING_GLB_CUT+") && " + PASSING_TMI_CUT;

PASSING_VBTFLIKE_CUTS = ("numberOfMatches > 1 && " # tracker muon, 2 stations
                         "muonID('GlobalMuonPromptTight') && "+ # and global, valid muon hits in global track, chi2/ndf < 10
                         "track.numberOfValidHits > 10 && track.hitPattern.numberOfValidPixelHits > 0 && "+ # track cuts
                         "abs(dB) < 0.2") # dxy w.r.t. PV or BS, 2mm

PASS_HLT_Mu3 = "(!triggerObjectMatchesByFilter('hltSingleMu3L3Filtered3').empty())"
PASS_HLT_L1DoubleMuOpen = "(!triggerObjectMatchesByFilter('hltDoubleMuLevel1PathL1OpenFiltered').empty())"
PASS_HLT_Mu0_Track0_Jpsi_New = ("(!triggerObjectMatchesByCollection('hltL3MuonCandidates::HLT').empty() && "+
                               "  triggerObjectMatchesByCollection('hltL3MuonCandidates::HLT').at(0).hasFilterLabel('hltMu0TrackJpsiTrackMassFiltered'))")
PASS_HLT_Mu0_Track0_Jpsi_ReDigi    = PASS_HLT_Mu0_Track0_Jpsi_New.replace("::HLT","::REDIGI")
PASS_HLT_Mu0_Track0_Jpsi_ReDigi362 = PASS_HLT_Mu0_Track0_Jpsi_New.replace("::HLT","::REDIGI36X")
PASS_HLT_Mu0_Track0_Jpsi = "( %s || %s || %s )" % ( PASS_HLT_Mu0_Track0_Jpsi_New, PASS_HLT_Mu0_Track0_Jpsi_ReDigi, PASS_HLT_Mu0_Track0_Jpsi_ReDigi362 )
PASS_HLT_Mu3_Track0_Jpsi = PASS_HLT_Mu0_Track0_Jpsi.replace("Mu0Track","Mu3Track")
PASS_HLT_1MU = "( %s || %s )" % ( PASS_HLT_Mu3, PASS_HLT_Mu0_Track0_Jpsi)
PASS_HLT_2MU = "( %s || %s )" % ( PASS_HLT_Mu3, PASS_HLT_L1DoubleMuOpen)
TAG_CUTS_1MU = "isGlobalMuon && " + TRACK_CUTS +' && '+ PASS_HLT_1MU
TAG_CUTS_2MU = "isGlobalMuon && " + TRACK_CUTS +' && '+ PASS_HLT_2MU

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
tkTracks = cms.EDProducer("ConcreteChargedCandidateProducer", 
    src = cms.InputTag("betterTracks"),      
    particleType = cms.string("mu+"),
) 

##    ____  _                  _    _    _                  
##   / ___|| |_ __ _ _ __   __| |  / \  | | ___  _ __   ___ 
##   \___ \| __/ _` | '_ \ / _` | / _ \ | |/ _ \| '_ \ / _ \
##    ___) | || (_| | | | | (_| |/ ___ \| | (_) | | | |  __/
##   |____/ \__\__,_|_| |_|\__,_/_/   \_\_|\___/|_| |_|\___|
##                                 
staOneValidHit = cms.EDProducer("TrackSelector",
    src = cms.InputTag("standAloneMuons","UpdatedAtVtx"), 
    cut = cms.string("numberOfValidHits > 0"),
)
staTracks = cms.EDProducer("ConcreteChargedCandidateProducer", 
    src = cms.InputTag("standAloneMuons","UpdatedAtVtx"), 
    particleType = cms.string("mu+"),
)
staTracksValidHits = staTracks.clone(src = "staOneValidHit")
                         
##    _____               _        
##   |_   _| __ __ _  ___| | _____ 
##     | || '__/ _` |/ __| |/ / __|
##     | || | | (_| | (__|   <\__ \
##     |_||_|  \__,_|\___|_|\_\___/
##                                 
##   
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
        abseta = cms.string("abs(eta)"),
    ),
    # choice of what defines a 'passing' probe
    flags = cms.PSet(),
    ## variables relative to the _tag_
    tagVariables = cms.PSet(
        pt  = cms.string("pt"),
        eta = cms.string("eta"),
    ),
    tagFlags = cms.PSet(
        L1MuOpen = cms.string("!triggerObjectMatchesByFilter('hltL1MuOpenL1Filtered0').empty()"),
        L2Mu0    = cms.string("!triggerObjectMatchesByFilter('hltL2Mu0L2Filtered0').empty()"),
        L2Mu3    = cms.string("!triggerObjectMatchesByFilter('hltSingleMu3L2Filtered3').empty()"),
        Mu3      = cms.string("!triggerObjectMatchesByFilter('hltSingleMu3L3Filtered3').empty()"),
        Mu5      = cms.string("!triggerObjectMatchesByFilter('hltSingleMu5L3Filtered5').empty()"),
        L1DoubleMuOpen  = cms.string("!triggerObjectMatchesByFilter('hltDoubleMuLevel1PathL1OpenFiltered').empty()"),
        Mu0_L1MuOpen    = cms.string("!triggerObjectMatchesByFilter('hltMu0L1MuOpenL3Filtered0').empty()"),
        Mu3_L1MuOpen    = cms.string("!triggerObjectMatchesByFilter('hltMu3L1MuOpenL3Filtered3').empty()"),
        Mu0_Track0_JPsi = cms.string(PASS_HLT_Mu0_Track0_Jpsi),
        Mu3_Track0_JPsi = cms.string(PASS_HLT_Mu3_Track0_Jpsi),
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
    staTracks   +
    staOneValidHit * staTracksValidHits +
    tagMuons1Mu + 
    tagMuons2Mu + 
    muMcMatch 
)

def addDiMuonSeparationVariables(process, sequence, treeProducer):
    from MuonAnalysis.TagAndProbe.nearbyMuonsInfo_cfi import nearbyMuonsInfo;
    tpp = treeProducer.tagProbePairs.moduleLabel
    if not hasattr(process, tpp+"DiMuonInfo"):
        setattr(process, tpp+"DiMuonInfo", nearbyMuonsInfo.clone(src = treeProducer.tagProbePairs))
        sequence.replace(getattr(process,tpp), getattr(process,tpp) + getattr(process, tpp+"DiMuonInfo"))
    if not hasattr(treeProducer, 'pairVariables'):
        treeProducer.pairVariables = cms.PSet()
        treeProducer.pairFlags     = cms.PSet()
    treeProducer.pairVariables.dphiVtxTimesQ = cms.InputTag(tpp+"DiMuonInfo", "dphiVtxTimesQ")
    treeProducer.pairVariables.drM2          = cms.InputTag(tpp+"DiMuonInfo", "drM2")
    treeProducer.pairVariables.dphiM2        = cms.InputTag(tpp+"DiMuonInfo", "dphiM2")
    treeProducer.pairVariables.distM2        = cms.InputTag(tpp+"DiMuonInfo", "distM2")
    treeProducer.pairVariables.drStaIn       = cms.InputTag(tpp+"DiMuonInfo", "drStaIn")
    treeProducer.pairVariables.dphiStaIn     = cms.InputTag(tpp+"DiMuonInfo", "dphiStaIn")

def allTPTreeProducers(process):
    for K in process.analyzers_().keys():
        V = getattr(process,K) # ; print K, V.type_()
        if V.type_() == "TagProbeFitTreeProducer":
            yield (K,V)


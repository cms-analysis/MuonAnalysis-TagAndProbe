import FWCore.ParameterSet.Config as cms

process = cms.Process("TagProbe")

process.load('FWCore.MessageService.MessageLogger_cfi')
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.source = cms.Source("PoolSource", 
    fileNames = cms.untracked.vstring(
        '/store/data/Run2010B/MuOnia/RECO/Nov4ReReco_v1/0179/FA21CF01-3BF1-DF11-A9A1-00215E22053A.root',
        '/store/data/Run2010B/MuOnia/RECO/Nov4ReReco_v1/0179/746FE647-2AF1-DF11-8D41-00215E21D56A.root',
        '/store/data/Run2010B/MuOnia/RECO/Nov4ReReco_v1/0179/6AC847A6-2FF1-DF11-9D58-00215E220F78.root',
        '/store/data/Run2010B/MuOnia/RECO/Nov4ReReco_v1/0179/665E9DFE-3AF1-DF11-8FB3-00215E22237C.root',
        '/store/data/Run2010B/MuOnia/RECO/Nov4ReReco_v1/0179/647A53EE-39F1-DF11-9BAA-00215E21D8E2.root',
        '/store/data/Run2010B/MuOnia/RECO/Nov4ReReco_v1/0179/46632181-31F1-DF11-AB04-00215E21DAB0.root',
        '/store/data/Run2010B/MuOnia/RECO/Nov4ReReco_v1/0179/32EBDC2C-29F1-DF11-A0CE-00215E93EFB8.root',
        '/store/data/Run2010B/MuOnia/RECO/Nov4ReReco_v1/0179/285421F0-2BF1-DF11-8688-00215E21DAB0.root',
        '/store/data/Run2010B/MuOnia/RECO/Nov4ReReco_v1/0178/FADEA9B4-18F1-DF11-82CD-00215E2211AC.root',
        '/store/data/Run2010B/MuOnia/RECO/Nov4ReReco_v1/0178/96DE97D2-16F1-DF11-A850-00215E21D462.root',
        '/store/data/Run2010B/MuOnia/RECO/Nov4ReReco_v1/0178/8457837A-17F1-DF11-929E-00215E21D972.root',
        '/store/data/Run2010B/MuOnia/RECO/Nov4ReReco_v1/0178/84454AB5-18F1-DF11-8DC1-00215E221248.root',
        '/store/data/Run2010B/MuOnia/RECO/Nov4ReReco_v1/0178/68D27AC3-18F1-DF11-928D-00215E221FBC.root',
        '/store/data/Run2010B/MuOnia/RECO/Nov4ReReco_v1/0178/0ADC2DCB-16F1-DF11-8C38-00215E21DAAA.root',
        '/store/data/Run2010B/MuOnia/RECO/Nov4ReReco_v1/0177/DEE0EE56-FBF0-DF11-8E70-00215E21DA98.root',
        '/store/data/Run2010B/MuOnia/RECO/Nov4ReReco_v1/0177/B8B651F1-F8F0-DF11-B78E-00215E221B48.root',
        '/store/data/Run2010B/MuOnia/RECO/Nov4ReReco_v1/0177/928089F6-F9F0-DF11-9429-00215E221506.root',
        '/store/data/Run2010B/MuOnia/RECO/Nov4ReReco_v1/0177/8E498B51-FBF0-DF11-A542-00215E21DD0E.root',
        '/store/data/Run2010B/MuOnia/RECO/Nov4ReReco_v1/0177/28A8C477-FDF0-DF11-82B9-00215E221FF2.root',
        '/store/data/Run2010B/MuOnia/RECO/Nov4ReReco_v1/0177/22CA2540-08F1-DF11-94D8-E41F13180DC8.root',
        '/store/data/Run2010B/MuOnia/RECO/Nov4ReReco_v1/0177/1EA26A07-0AF1-DF11-9A8B-00215E221B48.root',
        '/store/data/Run2010B/MuOnia/RECO/Nov4ReReco_v1/0177/12D8F056-FBF0-DF11-8595-00215E21DA98.root',
        '/store/data/Run2010B/MuOnia/RECO/Nov4ReReco_v1/0177/0632E104-0AF1-DF11-9571-00215E21D9F6.root',
        '/store/data/Run2010B/MuOnia/RECO/Nov4ReReco_v1/0177/00716ADD-0AF1-DF11-AD13-00215E22053A.root',
        '/store/data/Run2010B/MuOnia/RECO/Nov4ReReco_v1/0176/BA7734F2-EDF0-DF11-84B3-00215E21DB4C.root',
        '/store/data/Run2010B/MuOnia/RECO/Nov4ReReco_v1/0176/6A75AA1C-F0F0-DF11-A4CD-00215E22237C.root',
    ),
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )    


process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.GlobalTag.globaltag = cms.string('GR_R_38X_V14::All')

## ==== Fast Filters ====
process.goodVertexFilter = cms.EDFilter("VertexSelector",
    src = cms.InputTag("offlinePrimaryVertices"),
    cut = cms.string("!isFake && ndof > 4 && abs(z) <= 25 && position.Rho <= 2"),
    filter = cms.bool(True),
)
process.noScraping = cms.EDFilter("FilterOutScraping",
    applyfilter = cms.untracked.bool(True),
    debugOn = cms.untracked.bool(False), ## Or 'True' to get some per-event info
    numtrack = cms.untracked.uint32(10),
    thresh = cms.untracked.double(0.25)
)
process.fastFilter = cms.Sequence(process.goodVertexFilter + process.noScraping)

process.load("HLTrigger.HLTfilters.triggerResultsFilter_cfi")
process.triggerResultsFilter.triggerConditions = cms.vstring( 'HLT_Jet*' )
process.triggerResultsFilter.l1tResults = ''
process.triggerResultsFilter.throw = True
process.triggerResultsFilter.hltResults = cms.InputTag( "TriggerResults", "", "HLT" )
process.HLTMu   = process.triggerResultsFilter.clone(triggerConditions = [ 'HLT_Mu5_L2Mu0' ])
process.HLTBoth = process.triggerResultsFilter.clone(triggerConditions = [ 'HLT_Mu5_L2Mu0', 'HLT_Mu3_Track*_Jpsi*', 'HLT_Mu5_Track*_Jpsi*' ])

## ==== Merge CaloMuons and Tracks into the collection of reco::Muons  ====
from RecoMuon.MuonIdentification.calomuons_cfi import calomuons;
process.mergedMuons = cms.EDProducer("CaloMuonMerger",
    mergeTracks = cms.bool(True),
    muons     = cms.InputTag("muons"), 
    caloMuons = cms.InputTag("calomuons"),
    tracks    = cms.InputTag("generalTracks"),
    minCaloCompatibility = calomuons.minCaloCompatibility,
    ## Apply some minimal pt cut
    muonsCut     = cms.string("track.isNonnull"),
    caloMuonsCut = cms.string(""),
    tracksCut    = cms.string("pt > 2 || (abs(eta) > 1 && p > 2)"),
)

## ==== Trigger matching
process.load("MuonAnalysis.MuonAssociators.patMuonsWithTrigger_cff")
## with some customization
from MuonAnalysis.MuonAssociators.patMuonsWithTrigger_cff import *
changeRecoMuonInput(process, "mergedMuons")
#useL1MatchingWindowForSinglets(process) ## No longer used

from MuonAnalysis.TagAndProbe.common_variables_cff import *
process.load("MuonAnalysis.TagAndProbe.common_modules_cff")

process.tagMuons = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("patMuonsWithTrigger"),
    cut = cms.string("isGlobalMuon  && pt > 3 && !triggerObjectMatchesByCollection('hltL3MuonCandidates').empty()"),
)

process.oneTag  = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("tagMuons"), minNumber = cms.uint32(1))

process.probeMuons = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("patMuonsWithTrigger"),
    cut = cms.string("track.isNonnull && (!triggerObjectMatchesByCollection('hltMuTrackJpsiCtfTrackCands').empty() || !triggerObjectMatchesByCollection('hltL2MuonCandidates').empty())"),
)

process.tpPairs = cms.EDProducer("CandViewShallowCloneCombiner",
    cut = cms.string('2.8 < mass < 3.4'),
    decay = cms.string('tagMuons@+ probeMuons@-')
)
process.onePair = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("tpPairs"), minNumber = cms.uint32(1))

process.tpTree = cms.EDAnalyzer("TagProbeFitTreeProducer",
    # choice of tag and probe pairs, and arbitration
    tagProbePairs = cms.InputTag("tpPairs"),
    arbitration   = cms.string("OneProbe"),
    # probe variables: all useful ones
    variables = cms.PSet(
        KinematicVariables, 
        MuonIDVariables, 
        TrackQualityVariables, 
        #L1Variables, L2Variables, L3Variables
    ),
    flags = cms.PSet(
       TrackQualityFlags,
       MuonIDFlags,
       HighPtTriggerFlags,
       LowPtTriggerFlagsPhysics,
       LowPtTriggerFlagsEfficienciesProbe,
       ## A few other flags
       Track_QTF  = cms.string("track.numberOfValidHits > 11 && track.hitPattern.pixelLayersWithMeasurement > 1 && track.normalizedChi2 < 4 && abs(dB) < 3 && abs(track.dz) < 30"),
       Track_VBTF = cms.string("track.numberOfValidHits > 10 && track.hitPattern.pixelLayersWithMeasurement > 0 && track.normalizedChi2 < 4 && abs(dB) < 0.2"),
       ## Acceptance definition
       Acc_JPsi = cms.string("(abs(eta) <= 1.3 && pt > 3.3) || (1.3 < abs(eta) <= 2.2 && p > 2.9) || (2.2 < abs(eta) <= 2.4  && pt > 0.8)"),
    ),
    tagVariables = cms.PSet(
        pt  = cms.string('pt'),
        eta = cms.string('eta'),
        nVertices = cms.InputTag("nverticesModule"),
    ),
    tagFlags     = cms.PSet(
        #LowPtTriggerFlagsPhysics,
        #LowPtTriggerFlagsEfficienciesTag,
    ),
    pairVariables = cms.PSet(
        dphiVtxTimesQ = cms.InputTag("tagProbeSeparation", "dphiVtxTimesQ"),
        drM2          = cms.InputTag("tagProbeSeparation", "drM2"),
        dphiM2        = cms.InputTag("tagProbeSeparation", "dphiM2"),
        distM2        = cms.InputTag("tagProbeSeparation", "distM2"),
        #drStaIn       = cms.InputTag("tagProbeSeparation", "drStaIn"),
        #dphiStaIn     = cms.InputTag("tagProbeSeparation", "dphiStaIn"),
    ),
    pairFlags = cms.PSet(),
    isMC           = cms.bool(False),
    addRunLumiInfo = cms.bool(True),
)
process.tnpSimpleSequence = cms.Sequence(
    process.tagMuons +
    process.oneTag     +
    process.nverticesModule +
    process.probeMuons +
    process.tpPairs    +
    process.onePair    +
    process.tagProbeSeparation +
    process.tpTree
)

process.tagAndProbe = cms.Path( 
    process.fastFilter +
    process.HLTBoth    +
    process.mergedMuons                 *
    process.patMuonsWithTriggerSequence *
    process.tnpSimpleSequence
)

##    _____               _    _             
##   |_   _| __ __ _  ___| | _(_)_ __   __ _ 
##     | || '__/ _` |/ __| |/ / | '_ \ / _` |
##     | || | | (_| | (__|   <| | | | | (_| |
##     |_||_|  \__,_|\___|_|\_\_|_| |_|\__, |
##                                     |___/ 

## Then make another collection for standalone muons, using standalone track to define the 4-momentum
process.muonsSta = cms.EDProducer("RedefineMuonP4FromTrack",
    src   = cms.InputTag("muons"),
    track = cms.string("outer"),
)
## Match to trigger, to measure the efficiency of HLT tracking
from PhysicsTools.PatAlgos.tools.helpers import *
process.patMuonsWithTriggerSequenceSta = cloneProcessingSnippet(process, process.patMuonsWithTriggerSequence, "Sta")
process.muonMatchHLTL2Sta.maxDeltaR = 0.5
process.muonMatchHLTL3Sta.maxDeltaR = 0.5
massSearchReplaceAnyInputTag(process.patMuonsWithTriggerSequenceSta, "mergedMuons", "muonsSta")

## Define probes and T&P pairs
process.probeMuonsSta = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("patMuonsWithTriggerSta"),
    cut = cms.string("outerTrack.isNonnull && !triggerObjectMatchesByCollection('hltL2MuonCandidates').empty()"), 
)
process.tpPairsSta = process.tpPairs.clone(decay = "tagMuons@+ probeMuonsSta@-", cut = "2 < mass < 5")

process.onePairSta = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("tpPairsSta"), minNumber = cms.uint32(1))

## Now I have to define the passing probes for tracking
process.tkTracks = cms.EDProducer("ConcreteChargedCandidateProducer", 
    src = cms.InputTag("generalTracks"),      
    particleType = cms.string("mu+"),
)
## Then filter out the J/Psi's, to compute fake matching rate
process.tkTracksNoJPsi = cms.EDProducer("CandidateResonanceInefficiencyCreator",
    src = cms.InputTag("tkTracks"),
    tags = cms.InputTag("tagMuons"),
    mass    = cms.double(3.096),
    massMin = cms.double(2.85), ## Should cut away
    massMax = cms.double(3.25), ## 99.5% of signal
    onlyBestMatch = cms.bool(True),
    outputMode = cms.string("RefToBaseVector"),
)

process.staToTkMatch = cms.EDProducer("MatcherUsingTracks",
    src     = cms.InputTag("probeMuonsSta"),
    matched = cms.InputTag("tkTracks"),  
    algorithm = cms.string("byDirectComparison"), 
    srcTrack     = cms.string("muon"),    srcState = cms.string("atVertex"), 
    matchedTrack = cms.string("tracker"), matchedState = cms.string("atVertex"),
    maxDeltaR        = cms.double(1.),   # large range in DR (we can tighten it later)
    maxDeltaEta      = cms.double(0.4),  # small in eta, which is more precise
    maxDeltaLocalPos = cms.double(100),
    maxDeltaPtRel    = cms.double(5),   # |pt(sta) - pt(tk)|/pt(tk)
    sortBy           = cms.string("deltaR"),
)
process.staToTkMatchNoJPsi = process.staToTkMatch.clone(matched = 'tkTracksNoJPsi')
process.staPassingTk = cms.EDProducer("MatchedCandidateSelector",
    src   = cms.InputTag("probeMuonsSta"),
    match = cms.InputTag("staToTkMatch"),
)
process.staPassingTkNoJPsi = process.staPassingTk.clone(match = 'staToTkMatchNoJPsi')

process.tpTreeSta = process.tpTree.clone(
    tagProbePairs = "tpPairsSta",
    variables = cms.PSet(
        KinematicVariables, 
        #TriggerVariables, 
        ## extra standalone muon quality variables
        outerHits      = cms.string("outerTrack.hitPattern.numberOfHits"),
        outerValidHits = cms.string("outerTrack.numberOfValidHits"),
        outerStationsAny   = cms.string("outerTrack.hitPattern.muonStationsWithAnyHits"),
        outerStationsValid = cms.string("outerTrack.hitPattern.muonStationsWithValidHits"),
        ## track matching variables
        tk_deltaR     = cms.InputTag("staToTkMatch","deltaR"),
        tk_deltaPtRel = cms.InputTag("staToTkMatch","deltaPtRel"),
        tk_deltaEta   = cms.InputTag("staToTkMatch","deltaEta"),
        tk_deltaPhi   = cms.InputTag("staToTkMatch","deltaPhi"),
        tk_deltaR_NoJPsi     = cms.InputTag("staToTkMatchNoJPsi","deltaR"),
        tk_deltaPtRel_NoJPsi = cms.InputTag("staToTkMatchNoJPsi","deltaPtRel"),
        tk_deltaEta_NoJPsi   = cms.InputTag("staToTkMatchNoJPsi","deltaEta"),
        tk_deltaPhi_NoJPsi   = cms.InputTag("staToTkMatchNoJPsi","deltaPhi"),
    ),
    flags = cms.PSet(
        #LowPtTriggerFlagsPhysics,
        #LowPtTriggerFlagsEfficienciesProbe,
        hasTrack       = cms.InputTag("staPassingTk"),
        hasTrackNoJPsi = cms.InputTag("staPassingTkNoJPsi"),
    ),
    pairVariables = cms.PSet(),
    pairFlags     = cms.PSet(),
)
#process.tpTreeSta.variables.l1pt = process.tpTreeSta.variables.l1pt.value().replace("muonL1Info","muonL1InfoSta")
#process.tpTreeSta.variables.l1q  = process.tpTreeSta.variables.l1q.value( ).replace("muonL1Info","muonL1InfoSta")
#process.tpTreeSta.variables.l1dr = process.tpTreeSta.variables.l1dr.value().replace("muonL1Info","muonL1InfoSta")

process.tnpSimpleSequenceSta = cms.Sequence(
    process.tagMuons +
    process.oneTag     +
    process.nverticesModule +
    process.probeMuonsSta +
    ( process.tkTracks       * process.staToTkMatch       * process.staPassingTk       +
      process.tkTracksNoJPsi * process.staToTkMatchNoJPsi * process.staPassingTkNoJPsi ) +
    process.tpPairsSta      +
    process.onePairSta      +
    process.tpTreeSta
)

process.tagAndProbeSta = cms.Path( 
    process.fastFilter +
    process.HLTMu      +
    process.muonsSta                       +
    process.patMuonsWithTriggerSequenceSta +
    process.tnpSimpleSequenceSta
)

process.TFileService = cms.Service("TFileService", fileName = cms.string("tnpJPsi_Data.root"))

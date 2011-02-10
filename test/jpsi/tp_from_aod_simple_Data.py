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
process.GlobalTag.globaltag = cms.string('GR_R_39X_V5::All')

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
process.HLTMu   = process.triggerResultsFilter.clone(triggerConditions = [ 'HLT_Mu*_L2Mu0' ])
process.HLTBoth = process.triggerResultsFilter.clone(triggerConditions = [ 'HLT_Mu*_L2Mu0', 'HLT_Mu3_Track*_Jpsi*', 'HLT_Mu5_Track*_Jpsi*' ])

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

import MuonAnalysis.TagAndProbe.expectedHitsComputer_cfi
process.expectedHitsMu = MuonAnalysis.TagAndProbe.expectedHitsComputer_cfi.expectedHitsComputer.clone()
process.expectedHitsMu.inputColl  = cms.InputTag("mergedMuons")


## ==== Trigger matching
process.load("MuonAnalysis.MuonAssociators.patMuonsWithTrigger_cff")
## with some customization
from MuonAnalysis.MuonAssociators.patMuonsWithTrigger_cff import *
changeRecoMuonInput(process, "mergedMuons")
#useL1MatchingWindowForSinglets(process) ## No longer used
process.patMuonsWithoutTrigger.userData.userInts.src = cms.VInputTag(
    cms.InputTag('expectedHitsMu','in'),
    cms.InputTag('expectedHitsMu','out')
)


from MuonAnalysis.TagAndProbe.common_variables_cff import *
process.load("MuonAnalysis.TagAndProbe.common_modules_cff")

process.tagMuons = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("patMuonsWithTrigger"),
    cut = cms.string("(isGlobalMuon || numberOfMatchedStations > 1 ) && pt > 3 && !triggerObjectMatchesByCollection('hltL3MuonCandidates').empty()"),
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
        AllVariables,
        dxyPVdzmin = cms.InputTag("muonDxyPVdzmin","dxyPVdzmin"),
    ),
    flags = cms.PSet(
       TrackQualityFlags,
       MuonIDFlags,
       HighPtTriggerFlags,
       LowPtTriggerFlagsPhysics,
       LowPtTriggerFlagsEfficienciesProbe,
       ## ParticleFlow
       PF = cms.InputTag("muonsPassingPF"),
       ## A few other flags
       Track_QTF  = cms.string("track.numberOfValidHits > 11 && track.hitPattern.pixelLayersWithMeasurement > 1 && track.normalizedChi2 < 4 && abs(dB) < 3 && abs(track.dz) < 30"),
       Track_VBTF = cms.string("track.numberOfValidHits > 10 && track.hitPattern.pixelLayersWithMeasurement > 0 && abs(dB) < 0.2"),
       ## Acceptance definition
       Acc_JPsi = cms.string("(abs(eta) <= 1.3 && pt > 3.3) || (1.3 < abs(eta) <= 2.2 && p > 2.9) || (2.2 < abs(eta) <= 2.4  && pt > 0.8)"),
    ),
    tagVariables = cms.PSet(
        pt  = cms.string('pt'),
        eta = cms.string('eta'),
        nVertices = cms.InputTag("nverticesModule"),
    ),
    tagFlags     = cms.PSet(
        LowPtTriggerFlagsPhysics,
        LowPtTriggerFlagsEfficienciesTag,
    ),
    pairVariables = cms.PSet(
        pt = cms.string("pt"),
        dphiVtxTimesQ = cms.InputTag("tagProbeSeparation", "dphiVtxTimesQ"),
        drM1          = cms.InputTag("tagProbeSeparation", "drM1"),
        dphiM1        = cms.InputTag("tagProbeSeparation", "dphiM1"),
        distM1        = cms.InputTag("tagProbeSeparation", "distM1"),
        drM2          = cms.InputTag("tagProbeSeparation", "drM2"),
        dphiM2        = cms.InputTag("tagProbeSeparation", "dphiM2"),
        distM2        = cms.InputTag("tagProbeSeparation", "distM2"),
        drVtx         = cms.InputTag("tagProbeSeparation", "drVtx"),
    ),
    pairFlags = cms.PSet(),
    isMC           = cms.bool(False),
    addRunLumiInfo = cms.bool(True),
)
process.tnpSimpleSequence = cms.Sequence(
    process.tagMuons +
    process.oneTag     +
    process.probeMuons +
    process.tpPairs    +
    process.onePair    +
    process.muonDxyPVdzmin +
    process.nverticesModule +
    process.tagProbeSeparation +
    process.muonsPassingPF +
    process.tpTree
)

process.tagAndProbe = cms.Path( 
    process.fastFilter +
    process.HLTBoth    +
    process.mergedMuons                 *
    process.expectedHitsMu              *
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

process.tpTreeSta = process.tpTree.clone(
    tagProbePairs = "tpPairsSta",
    variables = cms.PSet(
        KinematicVariables, 
        ## track matching variables
        tk_deltaR     = cms.InputTag("staToTkMatch","deltaR"),
        tk_deltaEta   = cms.InputTag("staToTkMatch","deltaEta"),
        tk_deltaR_NoJPsi     = cms.InputTag("staToTkMatchNoJPsi","deltaR"),
        tk_deltaEta_NoJPsi   = cms.InputTag("staToTkMatchNoJPsi","deltaEta"),
        tk_deltaR_NoBestJPsi     = cms.InputTag("staToTkMatchNoBestJPsi","deltaR"),
        tk_deltaEta_NoBestJPsi   = cms.InputTag("staToTkMatchNoBestJPsi","deltaEta"),
    ),
    flags = cms.PSet(
        outerValidHits = cms.string("outerTrack.numberOfValidHits > 0"),
        MuX_L2Mu0_L2   = LowPtTriggerFlagsEfficienciesProbe.MuX_L2Mu0_L2,
        TM  = cms.string("isTrackerMuon"),
        Glb = cms.string("isGlobalMuon"),
    ),
    tagVariables = cms.PSet(
        nVertices = cms.InputTag("nverticesModule"),
    ),
    tagFlags = cms.PSet(
        Mu0_L2Mu0_MU = LowPtTriggerFlagsEfficienciesTag.Mu0_L2Mu0_MU,
        Mu3_L2Mu0_MU = LowPtTriggerFlagsEfficienciesTag.Mu3_L2Mu0_MU,
        Mu5_L2Mu0_MU = LowPtTriggerFlagsEfficienciesTag.Mu5_L2Mu0_MU,
    ),
    pairVariables = cms.PSet(),
    pairFlags     = cms.PSet(),
)

process.tnpSimpleSequenceSta = cms.Sequence(
    process.tagMuons +
    process.oneTag     +
    process.probeMuonsSta +
    process.tpPairsSta      +
    process.onePairSta      +
    process.nverticesModule +
    process.staToTkMatchSequenceJPsi +
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

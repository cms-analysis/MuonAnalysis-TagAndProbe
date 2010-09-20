import FWCore.ParameterSet.Config as cms

process = cms.Process("TagProbe")

process.load('FWCore.MessageService.MessageLogger_cfi')
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.source = cms.Source("PoolSource", 
    fileNames = cms.untracked.vstring(
	'/store/data/Run2010A/MuOnia/RECO/v4/000/140/401/FC6A5623-1993-DF11-B37D-0030487C7828.root',
	'/store/data/Run2010A/MuOnia/RECO/v4/000/140/401/F202B792-1093-DF11-889A-001617C3B70E.root',
	'/store/data/Run2010A/MuOnia/RECO/v4/000/140/401/D8023811-0D93-DF11-8ECA-0030487CAF5E.root',
	'/store/data/Run2010A/MuOnia/RECO/v4/000/140/401/CC367D91-1093-DF11-B545-001617C3B6CC.root',
	'/store/data/Run2010A/MuOnia/RECO/v4/000/140/401/C0C63A7A-1393-DF11-9C57-003048F11114.root',
	'/store/data/Run2010A/MuOnia/RECO/v4/000/140/401/88C873F1-1193-DF11-BBFD-0030487A18F2.root',
	'/store/data/Run2010A/MuOnia/RECO/v4/000/140/401/7E56E31E-1993-DF11-B4C5-000423D98950.root',
	'/store/data/Run2010A/MuOnia/RECO/v4/000/140/401/66C7120F-1493-DF11-B150-0030487CAEAC.root',
	'/store/data/Run2010A/MuOnia/RECO/v4/000/140/401/343F2303-0693-DF11-BBE3-003048F1BF68.root',
	'/store/data/Run2010A/MuOnia/RECO/v4/000/140/401/3041B3C8-1B93-DF11-9C44-003048F110BE.root',
	'/store/data/Run2010A/MuOnia/RECO/v4/000/140/401/22CEE4C7-1B93-DF11-B7C7-003048F1182E.root',
    ),
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )    


process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.GlobalTag.globaltag = cms.string('GR_R_38X_V8::All')

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
useL1MatchingWindowForSinglets(process)

from MuonAnalysis.TagAndProbe.common_variables_cff import *
process.load("MuonAnalysis.TagAndProbe.common_modules_cff")

process.tagMuons = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("patMuonsWithTrigger"),
    cut = cms.string("isGlobalMuon"),
)

process.probeMuons = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("patMuonsWithTrigger"),
    cut = cms.string("track.isNonnull"),
)

process.tpPairs = cms.EDProducer("CandViewShallowCloneCombiner",
    cut = cms.string('2.6 < mass < 3.6'),
    decay = cms.string('tagMuons@+ probeMuons@-')
)

process.tpTree = cms.EDAnalyzer("TagProbeFitTreeProducer",
    # choice of tag and probe pairs, and arbitration
    tagProbePairs = cms.InputTag("tpPairs"),
    arbitration   = cms.string("OneProbe"),
    # probe variables: all useful ones
    variables = cms.PSet(
        KinematicVariables, 
        MuonIDVariables, 
        TrackQualityVariables, 
        L1Variables, L2Variables, L3Variables
    ),
    flags = cms.PSet(
       TrackQualityFlags,
       MuonIDFlags,
       HighPtTriggerFlags,
       LowPtTriggerFlagsPhysics,
       LowPtTriggerFlagsEfficiencies,
        ## A few other flags
        Track_QTF  = cms.string("track.numberOfValidHits > 11 && track.hitPattern.pixelLayersWithMeasurement > 1 && track.normalizedChi2 < 4 && abs(dB) < 3 && abs(track.dz) < 30"),
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
        LowPtTriggerFlagsEfficiencies,
    ),
    isMC           = cms.bool(False),
    addRunLumiInfo = cms.bool(True),
)

triggers4Mu = ['Mu3','Mu0_Track0_Jpsi','Mu3_Track0_Jpsi','Mu5_Track0_Jpsi','Mu3_Track3_Jpsi']
triggers4Tk = ['Mu3','L1DoubleMuOpen','L1DoubleMuOpen_Tight','L2DoubleMu0','Mu5_L2Mu0']
process.tpTree.tagFlags.unbias4Mu = cms.string(mkUnbiasCut(process.tpTree.tagFlags, triggers4Mu))
process.tpTree.tagFlags.unbias4Tk = cms.string(mkUnbiasCut(process.tpTree.tagFlags, triggers4Tk))

process.tnpSimpleSequence = cms.Sequence(
    process.tagMuons +
    process.nverticesModule +
    process.probeMuons +
    process.tpPairs    +
    process.tpTree
)


process.tagAndProbe = cms.Path( 
    process.fastFilter +
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
    cut = cms.string("outerTrack.isNonnull"), # no real cut now
)
process.tpPairsSta = process.tpPairs.clone(decay = "tagMuons@+ probeMuonsSta@-", cut = "2 < mass < 5")

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
        TriggerVariables, 
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
        LowPtTriggerFlagsPhysics,
        LowPtTriggerFlagsEfficiencies,
        hasTrack       = cms.InputTag("staPassingTk"),
        hasTrackNoJPsi = cms.InputTag("staPassingTkNoJPsi"),
    ),
)
process.tpTreeSta.variables.l1pt = process.tpTreeSta.variables.l1pt.value().replace("muonL1Info","muonL1InfoSta")
process.tpTreeSta.variables.l1q  = process.tpTreeSta.variables.l1q.value( ).replace("muonL1Info","muonL1InfoSta")
process.tpTreeSta.variables.l1dr = process.tpTreeSta.variables.l1dr.value().replace("muonL1Info","muonL1InfoSta")

process.tnpSimpleSequenceSta = cms.Sequence(
    process.tagMuons +
    process.nverticesModule +
    process.probeMuonsSta +
    ( process.tkTracks       * process.staToTkMatch       * process.staPassingTk       +
      process.tkTracksNoJPsi * process.staToTkMatchNoJPsi * process.staPassingTkNoJPsi ) +
    process.tpPairsSta      +
    process.tpTreeSta
)

process.tagAndProbeSta = cms.Path( 
    process.fastFilter +
    process.muonsSta                       +
    process.patMuonsWithTriggerSequenceSta +
    process.tnpSimpleSequenceSta
)

process.TFileService = cms.Service("TFileService", fileName = cms.string("tnpJPsi_Data.root"))

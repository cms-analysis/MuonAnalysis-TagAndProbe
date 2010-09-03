import FWCore.ParameterSet.Config as cms

process = cms.Process("TagProbe")

process.load('FWCore.MessageService.MessageLogger_cfi')
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.source = cms.Source("PoolSource", 
    fileNames = cms.untracked.vstring(
        '/store/data/Run2010A/Mu/RECO/v4/000/144/112/FE61BD50-CAB3-DF11-B7FD-001D09F29321.root',
        '/store/data/Run2010A/Mu/RECO/v4/000/144/112/FCEBD143-D1B3-DF11-854A-001D09F2423B.root',
        '/store/data/Run2010A/Mu/RECO/v4/000/144/112/FC5756A5-E3B3-DF11-8A7A-001617E30CE8.root',
        '/store/data/Run2010A/Mu/RECO/v4/000/144/112/F46E9FD9-01B4-DF11-BC69-001D09F2424A.root',
        '/store/data/Run2010A/Mu/RECO/v4/000/144/112/F0EFC873-E4B3-DF11-9E2F-0030487CD162.root',
        '/store/data/Run2010A/Mu/RECO/v4/000/144/112/EE66B7B6-0BB4-DF11-BC6C-001D09F2B30B.root',
        '/store/data/Run2010A/Mu/RECO/v4/000/144/112/E6F6D3D9-CBB3-DF11-8A0E-001D09F2AD4D.root',
        '/store/data/Run2010A/Mu/RECO/v4/000/144/112/E6F21B05-07B4-DF11-A21F-001D09F2447F.root',
        '/store/data/Run2010A/Mu/RECO/v4/000/144/112/E655B00D-C6B3-DF11-809B-001D09F2A49C.root',
        '/store/data/Run2010A/Mu/RECO/v4/000/144/112/E4DFA7CC-04B4-DF11-8C25-0030487CD6DA.root',
        '/store/data/Run2010A/Mu/RECO/v4/000/144/112/E0EEEDC8-F8B3-DF11-BBB7-0030487CD906.root',
        '/store/data/Run2010A/Mu/RECO/v4/000/144/112/D4BCA41B-ECB3-DF11-88F2-001D09F29849.root',
    ),
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )    

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.GlobalTag.globaltag = cms.string('GR_R_38X_V8::All')

##    __  __                       
##   |  \/  |_   _  ___  _ __  ___ 
##   | |\/| | | | |/ _ \| '_ \/ __|
##   | |  | | |_| | (_) | | | \__ \
##   |_|  |_|\__,_|\___/|_| |_|___/
##                                 
## ==== Merge CaloMuons and Tracks into the collection of reco::Muons  ====
from RecoMuon.MuonIdentification.calomuons_cfi import calomuons;
process.mergedMuons = cms.EDProducer("CaloMuonMerger",
    mergeTracks = cms.bool(True),
    muons     = cms.InputTag("muons"), 
    caloMuons = cms.InputTag("calomuons"),
    tracks    = cms.InputTag("generalTracks"),
    minCaloCompatibility = calomuons.minCaloCompatibility,
    ## Apply some minimal pt cut
    muonsCut     = cms.string("pt > 5 && track.isNonnull"),
    caloMuonsCut = cms.string("pt > 5"),
    tracksCut    = cms.string("pt > 5"),
)

## ==== Trigger matching
process.load("MuonAnalysis.MuonAssociators.patMuonsWithTrigger_cff")
## with some customization
process.muonMatchHLTL2.maxDeltaR = 0.5
process.muonMatchHLTL3.maxDeltaR = 0.1
from MuonAnalysis.MuonAssociators.patMuonsWithTrigger_cff import *
changeRecoMuonInput(process, "mergedMuons")

from MuonAnalysis.TagAndProbe.common_variables_cff import *
process.load("MuonAnalysis.TagAndProbe.common_modules_cff")

process.tagMuons = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("patMuonsWithTrigger"),
    cut = cms.string("pt > 15 && "+MuonIDFlags.VBTF.value()+" && !triggerObjectMatchesByFilter('hltSingleMu9L3Filtered9').empty()"),
)

process.probeMuons = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("patMuonsWithTrigger"),
    cut = cms.string("track.isNonnull"),  # no real cut now
)

process.tpPairs = cms.EDProducer("CandViewShallowCloneCombiner",
    cut = cms.string('40 < mass < 140'),
    decay = cms.string('tagMuons@+ probeMuons@-')
)

process.tpTree = cms.EDAnalyzer("TagProbeFitTreeProducer",
    # choice of tag and probe pairs, and arbitration
    tagProbePairs = cms.InputTag("tpPairs"),
    arbitration   = cms.string("OneProbe"),
    # probe variables: all useful ones
    variables = AllVariables,
    flags = cms.PSet(
       TrackQualityFlags,
       MuonIDFlags,
       HighPtTriggerFlags,
       ## Isolation
       Isol    = cms.string("(isolationR03.emEt + isolationR03.hadEt + isolationR03.sumPt)/pt < 0.15"), 
       IsolTk3 = cms.string("isolationR03.sumPt < 3"), 
    ),
    tagVariables = cms.PSet(
        nVertices = cms.InputTag("nverticesModule"),
    ),
    tagFlags = cms.PSet(),
    isMC           = cms.bool(False),
    addRunLumiInfo = cms.bool(True),
)

process.tnpSimpleSequence = cms.Sequence(
    process.tagMuons   *
    process.nverticesModule +
    process.probeMuons +
    process.tpPairs    *
    process.tpTree
)

from HLTrigger.HLTfilters.hltHighLevelDev_cfi import hltHighLevelDev
process.fastFilter = hltHighLevelDev.clone(HLTPaths = ['HLT_Mu9'], HLTPathsPrescales = [1])

process.tagAndProbe = cms.Path( 
    process.fastFilter                  +
    process.mergedMuons                 +
    process.patMuonsWithTriggerSequence +
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
process.tpPairsSta = process.tpPairs.clone(decay = "tagMuons@+ probeMuonsSta@-")

## Now I have to define the passing probes for tracking
process.tkTracks = cms.EDProducer("ConcreteChargedCandidateProducer", 
    src = cms.InputTag("generalTracks"),      
    particleType = cms.string("mu+"),
) 
process.staToTkMatch = cms.EDProducer("MatcherUsingTracks",
    src     = cms.InputTag("probeMuonsSta"),
    matched = cms.InputTag("tkTracks"),  
    algorithm = cms.string("byDirectComparison"), 
    srcTrack     = cms.string("muon"),    srcState = cms.string("atVertex"), 
    matchedTrack = cms.string("tracker"), matchedState = cms.string("atVertex"),
    maxDeltaR        = cms.double(0.3), 
    maxDeltaPtRel    = cms.double(2),   # |pt(sta) - pt(tk)|/pt(tk)
    maxDeltaLocalPos = cms.double(100),
    sortBy           = cms.string("deltaR"),
)
process.staPassingTk = cms.EDProducer("MatchedCandidateSelector",
    src   = cms.InputTag("probeMuonsSta"),
    match = cms.InputTag("staToTkMatch"),
)

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
    ),
    flags = cms.PSet(
        HighPtTriggerFlags, 
        hasTrack = cms.InputTag("staPassingTk"),
        L1DoubleMuOpen       = LowPtTriggerFlagsPhysics.L1DoubleMuOpen,
        L1DoubleMuOpen_Tight = LowPtTriggerFlagsPhysics.L1DoubleMuOpen_Tight,
        L2DoubleMu0          = LowPtTriggerFlagsPhysics.L2DoubleMu0,
    )
)
process.tpTreeSta.variables.l1pt = process.tpTreeSta.variables.l1pt.value().replace("muonL1Info","muonL1InfoSta")
process.tpTreeSta.variables.l1q  = process.tpTreeSta.variables.l1q.value( ).replace("muonL1Info","muonL1InfoSta")
process.tpTreeSta.variables.l1dr = process.tpTreeSta.variables.l1dr.value().replace("muonL1Info","muonL1InfoSta")
process.tpTreeSta.tagFlags = process.tpTreeSta.flags.clone(hasTrack = cms.string(""))

process.tnpSimpleSequenceSta = cms.Sequence(
    process.tagMuons        +
    process.nverticesModule +
    process.probeMuonsSta   +
    ( process.tkTracks * process.staToTkMatch * process.staPassingTk ) +
    process.tpPairsSta      +
    process.tpTreeSta
)

process.tagAndProbeSta = cms.Path( 
    process.fastFilter                     +
    process.muonsSta                       +
    process.patMuonsWithTriggerSequenceSta +
    process.tnpSimpleSequenceSta
)

process.TFileService = cms.Service("TFileService", fileName = cms.string("tnpZ.root"))

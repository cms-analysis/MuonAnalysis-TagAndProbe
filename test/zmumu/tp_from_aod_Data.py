import FWCore.ParameterSet.Config as cms

process = cms.Process("TagProbe")

process.load('FWCore.MessageService.MessageLogger_cfi')
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.source = cms.Source("PoolSource", 
    fileNames = cms.untracked.vstring(
	'/store/relval/CMSSW_3_9_5/Mu/RECO/GR_R_39X_V3_RelVal_wzMu2010B_realData-v1/0009/5ADD41D8-AEFA-DF11-8712-001A92971BD6.root',
	'/store/relval/CMSSW_3_9_5/Mu/RECO/GR_R_39X_V3_RelVal_wzMu2010B_realData-v1/0009/5618434E-B0FA-DF11-90B9-00261894389E.root',
	'/store/relval/CMSSW_3_9_5/Mu/RECO/GR_R_39X_V3_RelVal_wzMu2010B_realData-v1/0009/46336845-B1FA-DF11-A075-001A92971AA8.root',
	'/store/relval/CMSSW_3_9_5/Mu/RECO/GR_R_39X_V3_RelVal_wzMu2010B_realData-v1/0009/3E59C4E7-A8FA-DF11-BCE2-0018F3D09630.root',
	'/store/relval/CMSSW_3_9_5/Mu/RECO/GR_R_39X_V3_RelVal_wzMu2010B_realData-v1/0009/3629706C-A6FA-DF11-AABF-003048679044.root',
	'/store/relval/CMSSW_3_9_5/Mu/RECO/GR_R_39X_V3_RelVal_wzMu2010B_realData-v1/0009/0CE0A6C3-B3FA-DF11-AFEF-0018F3D0962A.root',
	'/store/relval/CMSSW_3_9_5/Mu/RECO/GR_R_39X_V3_RelVal_wzMu2010B_realData-v1/0009/0C8A93EF-A5FA-DF11-BCF7-0018F3D0961E.root',
	'/store/relval/CMSSW_3_9_5/Mu/RECO/GR_R_39X_V3_RelVal_wzMu2010B_realData-v1/0009/005823E0-A9FA-DF11-8949-002618943849.root',
	'/store/relval/CMSSW_3_9_5/Mu/RECO/GR_R_39X_V3_RelVal_wzMu2010B_realData-v1/0000/FE39D06D-99FA-DF11-9576-0018F3D096E4.root',
	'/store/relval/CMSSW_3_9_5/Mu/RECO/GR_R_39X_V3_RelVal_wzMu2010B_realData-v1/0000/F89C6E54-99FA-DF11-AF09-00304867BFC6.root',
	'/store/relval/CMSSW_3_9_5/Mu/RECO/GR_R_39X_V3_RelVal_wzMu2010B_realData-v1/0000/F6795867-99FA-DF11-9513-0018F3D09620.root',
	'/store/relval/CMSSW_3_9_5/Mu/RECO/GR_R_39X_V3_RelVal_wzMu2010B_realData-v1/0000/F6325427-9DFA-DF11-82A0-0018F3D096BA.root',
	'/store/relval/CMSSW_3_9_5/Mu/RECO/GR_R_39X_V3_RelVal_wzMu2010B_realData-v1/0000/F2E40049-99FA-DF11-99D9-00304867D838.root',
	'/store/relval/CMSSW_3_9_5/Mu/RECO/GR_R_39X_V3_RelVal_wzMu2010B_realData-v1/0000/EEE27061-99FA-DF11-9C06-001A92811724.root',
	'/store/relval/CMSSW_3_9_5/Mu/RECO/GR_R_39X_V3_RelVal_wzMu2010B_realData-v1/0000/EEE27061-99FA-DF11-84D9-001A92811724.root',
	'/store/relval/CMSSW_3_9_5/Mu/RECO/GR_R_39X_V3_RelVal_wzMu2010B_realData-v1/0000/E837E490-A7FA-DF11-B97C-0018F3D09604.root',
	#'/store/data/Run2010B/Mu/RECO/PromptReco-v2/000/149/181/F62B224C-09E4-DF11-8813-003048F117B6.root',
	#'/store/data/Run2010B/Mu/RECO/PromptReco-v2/000/149/181/E8E4D990-1BE4-DF11-9A0D-001D09F242EF.root',
	#'/store/data/Run2010B/Mu/RECO/PromptReco-v2/000/149/181/E25D7CF7-1CE4-DF11-AC3C-001D09F251B8.root',
	#'/store/data/Run2010B/Mu/RECO/PromptReco-v2/000/149/181/C45B1E12-FEE3-DF11-81BA-0030487C60AE.root',
	#'/store/data/Run2010B/Mu/RECO/PromptReco-v2/000/149/181/8840D649-EDE3-DF11-9429-003048D2C1C4.root',
	#'/store/data/Run2010B/Mu/RECO/PromptReco-v2/000/149/181/6AFEA573-EAE3-DF11-B597-0030487CD812.root',
	#'/store/data/Run2010B/Mu/RECO/PromptReco-v2/000/149/181/62DA1B87-14E4-DF11-B9CC-001D09F24F65.root',
	#'/store/data/Run2010B/Mu/RECO/PromptReco-v2/000/149/181/546FB2D5-E4E3-DF11-B43F-003048D2C0F0.root',
	#'/store/data/Run2010B/Mu/RECO/PromptReco-v2/000/149/181/46908DDD-E4E3-DF11-8E3A-003048D37580.root',
	#'/store/data/Run2010B/Mu/RECO/PromptReco-v2/000/149/181/4431849D-EEE3-DF11-925E-001D09F25217.root',
	#'/store/data/Run2010B/Mu/RECO/PromptReco-v2/000/149/181/3C95E446-EDE3-DF11-ADFD-003048D2C108.root',
	#'/store/data/Run2010B/Mu/RECO/PromptReco-v2/000/149/181/32671F77-19E4-DF11-A416-001D09F2462D.root',
	#'/store/data/Run2010B/Mu/RECO/PromptReco-v2/000/149/181/18682728-F2E3-DF11-B7F2-001617E30CD4.root',
	#'/store/data/Run2010B/Mu/RECO/PromptReco-v2/000/149/181/02D5E9A6-F5E3-DF11-9FCA-0030487CD6B4.root',
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
    #mergeCaloMuons = cms.bool(False), # AOD
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
    cut = cms.string("pt > 15 && "+MuonIDFlags.VBTF.value()+" && (!triggerObjectMatchesByPath('HLT_Mu9').empty() || !triggerObjectMatchesByPath('HLT_Mu15_v1').empty())"),
)

process.oneTag  = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("tagMuons"), minNumber = cms.uint32(1))

process.probeMuons = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("patMuonsWithTrigger"),
    cut = cms.string("track.isNonnull"),  # no real cut now
)

process.tpPairs = cms.EDProducer("CandViewShallowCloneCombiner",
    cut = cms.string('60 < mass < 140'),
    decay = cms.string('tagMuons@+ probeMuons@-')
)
process.onePair = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("tpPairs"), minNumber = cms.uint32(1))


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
       ## ParticleFlow
       PF = cms.InputTag("muonsPassingPF"),
       ## A few other flags
       Track_QTF  = cms.string("track.numberOfValidHits > 11 && track.hitPattern.pixelLayersWithMeasurement > 1 && track.normalizedChi2 < 4 && abs(dB) < 3 && abs(track.dz) < 30"),
       Track_VBTF = cms.string("track.numberOfValidHits > 10 && track.hitPattern.pixelLayersWithMeasurement > 0 && abs(dB) < 0.2"),
    ),
    tagVariables = cms.PSet(
        nVertices = cms.InputTag("nverticesModule"),
    ),
    tagFlags = cms.PSet(),
    pairVariables = cms.PSet(
        nJets15 = cms.InputTag("njets15Module"),
        nJets30 = cms.InputTag("njets30Module"),
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
    process.nverticesModule +
    process.njets15Module +
    process.njets30Module +
    process.muonsPassingPF +
    process.tpTree
)

process.tagAndProbe = cms.Path( 
    process.fastFilter +
    process.mergedMuons                 *
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

process.onePairSta = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("tpPairsSta"), minNumber = cms.uint32(1))

process.staToTkMatch.maxDeltaR     = 0.3
process.staToTkMatch.maxDeltaPtRel = 2.
process.staToTkMatchNoZ.maxDeltaR     = 0.3
process.staToTkMatchNoZ.maxDeltaPtRel = 2.

process.tpTreeSta = process.tpTree.clone(
    tagProbePairs = "tpPairsSta",
    variables = cms.PSet(
        KinematicVariables, 
        ## track matching variables
        tk_deltaR     = cms.InputTag("staToTkMatch","deltaR"),
        tk_deltaEta   = cms.InputTag("staToTkMatch","deltaEta"),
        tk_deltaR_NoZ   = cms.InputTag("staToTkMatchNoZ","deltaR"),
        tk_deltaEta_NoZ = cms.InputTag("staToTkMatchNoZ","deltaEta"),
    ),
    flags = cms.PSet(
        outerValidHits = cms.string("outerTrack.numberOfValidHits > 0"),
        TM  = cms.string("isTrackerMuon"),
        Glb = cms.string("isGlobalMuon"),
    ),
)
process.tpTreeSta.pairVariables.nJets15 = "njets15ModuleSta"
process.tpTreeSta.pairVariables.nJets30 = "njets30ModuleSta"
process.njets15ModuleSta = process.njets15Module.clone(pairs = "tpPairsSta")
process.njets30ModuleSta = process.njets30Module.clone(pairs = "tpPairsSta")

process.tnpSimpleSequenceSta = cms.Sequence(
    process.tagMuons +
    process.oneTag     +
    process.probeMuonsSta   +
    process.tpPairsSta      +
    process.onePairSta      +
    process.nverticesModule +
    process.staToTkMatchSequenceZ +
    process.njets15ModuleSta +
    process.njets30ModuleSta +
    process.tpTreeSta
)

process.tagAndProbeSta = cms.Path( 
    process.fastFilter +
    process.muonsSta                       +
    process.patMuonsWithTriggerSequenceSta +
    process.tnpSimpleSequenceSta
)

process.TFileService = cms.Service("TFileService", fileName = cms.string("tnpZ.root"))

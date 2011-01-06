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
process.GlobalTag.globaltag = cms.string('GR_R_39X_V4::All')

##    __  __                       
##   |  \/  |_   _  ___  _ __  ___ 
##   | |\/| | | | |/ _ \| '_ \/ __|
##   | |  | | |_| | (_) | | | \__ \
##   |_|  |_|\__,_|\___/|_| |_|___/
##                                 
## ==== Merge CaloMuons and Tracks into the collection of reco::Muons  ====
from RecoMuon.MuonIdentification.calomuons_cfi import calomuons;
process.mergedMuons = cms.EDProducer("CaloMuonMerger",
    mergeCaloMuons = cms.bool(False), # AOD
    mergeTracks    = cms.bool(True),
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

process.probeMuons = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("patMuonsWithTrigger"),
    cut = cms.string("track.isNonnull"),  # no real cut now
)

process.tpPairs = cms.EDProducer("CandViewShallowCloneCombiner",
    cut = cms.string('60 < mass < 140'),
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
    pairVariables = cms.PSet(
        dZ      = cms.string("daughter(0).vz - daughter(1).vz"),
        nJets15 = cms.InputTag("njets15Module"),
        nJets30 = cms.InputTag("njets30Module"),
    ),
    pairFlags = cms.PSet(),
    isMC           = cms.bool(False),
    addRunLumiInfo = cms.bool(True),
)

process.tnpSimpleSequence = cms.Sequence(
    process.tagMuons   *
    process.probeMuons +
    process.tpPairs    +
    process.nverticesModule +
    process.njets15Module +
    process.njets30Module +
    process.tpTree
)

process.tagAndProbe = cms.Path( 
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
process.staToTkNoZMatch = process.staToTkMatch.clone(matched = "tkTracksNoZ")
process.staPassingTk = cms.EDProducer("MatchedCandidateSelector",
    src   = cms.InputTag("probeMuonsSta"),
    match = cms.InputTag("staToTkMatch"),
)
process.staPassingTkNoZ = process.staPassingTk.clone(match = "staToTkNoZMatch")

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
        ## track matching variables for fake rates
        tkNoZ_deltaR     = cms.InputTag("staToTkNoZMatch","deltaR"),
        tkNoZ_deltaPtRel = cms.InputTag("staToTkNoZMatch","deltaPtRel"),
        tkNoZ_deltaEta   = cms.InputTag("staToTkNoZMatch","deltaEta"),
        tkNoZ_deltaPhi   = cms.InputTag("staToTkNoZMatch","deltaPhi"),
    ),
    flags = cms.PSet(
        #HighPtTriggerFlags, 
        Glb      = MuonIDFlags.Glb,
        TM       = MuonIDFlags.TM,
        TMA      = MuonIDFlags.TMA,
        hasTrack    = cms.InputTag("staPassingTk"),
        hasTrackNoZ = cms.InputTag("staPassingTkNoZ"),
        #L1DoubleMuOpen       = LowPtTriggerFlagsPhysics.L1DoubleMuOpen,
        #L1DoubleMuOpen_Tight = LowPtTriggerFlagsPhysics.L1DoubleMuOpen_Tight,
        #L2DoubleMu0          = LowPtTriggerFlagsPhysics.L2DoubleMu0,
    ),
)
process.tpTreeSta.variables.l1pt = process.tpTreeSta.variables.l1pt.value().replace("muonL1Info","muonL1InfoSta")
process.tpTreeSta.variables.l1q  = process.tpTreeSta.variables.l1q.value( ).replace("muonL1Info","muonL1InfoSta")
process.tpTreeSta.variables.l1dr = process.tpTreeSta.variables.l1dr.value().replace("muonL1Info","muonL1InfoSta")
process.tpTreeSta.tagFlags = process.tpTreeSta.flags.clone(hasTrack = cms.string(""), hasTrackNoZ  = cms.string(""))
process.tpTreeSta.pairVariables.nJets15 = "njets15ModuleSta"
process.tpTreeSta.pairVariables.nJets30 = "njets30ModuleSta"
process.njets15ModuleSta = process.njets15Module.clone(pairs = "tpPairsSta")
process.njets30ModuleSta = process.njets30Module.clone(pairs = "tpPairsSta")

process.tnpSimpleSequenceSta = cms.Sequence(
    process.tagMuons        +
    process.probeMuonsSta   +
    ( process.tkTracks    * process.staToTkMatch    * process.staPassingTk    ) +
    ( process.tkTracksNoZ * process.staToTkNoZMatch * process.staPassingTkNoZ ) +
    process.tpPairsSta      +
    process.nverticesModule +
    process.njets15ModuleSta +
    process.njets30ModuleSta +
    process.tpTreeSta
)

process.tagAndProbeSta = cms.Path( 
    process.muonsSta                       +
    process.patMuonsWithTriggerSequenceSta +
    process.tnpSimpleSequenceSta
)

process.TFileService = cms.Service("TFileService", fileName = cms.string("tnpZ.root"))

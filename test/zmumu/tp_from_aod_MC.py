import FWCore.ParameterSet.Config as cms

process = cms.Process("TagProbe")

process.load('FWCore.MessageService.MessageLogger_cfi')
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.source = cms.Source("PoolSource", 
    fileNames = cms.untracked.vstring(
	'/store/relval/CMSSW_5_2_3/RelValZMM/GEN-SIM-RECO/START52_V5-v1/0043/A29B9025-0E7A-E111-97E7-001A928116DE.root',
	'/store/relval/CMSSW_5_2_3/RelValZMM/GEN-SIM-RECO/START52_V5-v1/0043/5CAA0235-0F7A-E111-BA3E-0018F3D09690.root',
	'/store/relval/CMSSW_5_2_3/RelValZMM/GEN-SIM-RECO/START52_V5-v1/0043/1011EE9E-2B7A-E111-9349-0018F3D0970C.root',
	'/store/relval/CMSSW_5_2_3/RelValZMM/GEN-SIM-RECO/START52_V5-v1/0043/0E187509-0D7A-E111-8FA3-001A928116C2.root',
    ),
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(20000) )    

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.GlobalTag.globaltag = cms.string('START52_V5::All')

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
    mergeCaloMuons = cms.bool(False), # AOD
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
process.muonMatchHLTL2.maxDeltaR = 0.3 # Zoltan tuning - it was 0.5
process.muonMatchHLTL3.maxDeltaR = 0.1
from MuonAnalysis.MuonAssociators.patMuonsWithTrigger_cff import *
changeRecoMuonInput(process, "mergedMuons")
changeTriggerProcessName(process, "*") # auto-guess

from MuonAnalysis.TagAndProbe.common_variables_cff import *
process.load("MuonAnalysis.TagAndProbe.common_modules_cff")

process.tagMuons = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("patMuonsWithTrigger"),
    cut = cms.string("pt > 15 && "+MuonIDFlags.VBTF.value()+" && !triggerObjectMatchesByCollection('hltL3MuonCandidates').empty()"),
)

process.probeMuons = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("patMuonsWithTrigger"),
    cut = cms.string("track.isNonnull"),  # no real cut now
)

process.tpPairs = cms.EDProducer("CandViewShallowCloneCombiner",
    cut = cms.string('60 < mass < 140'),
    decay = cms.string('tagMuons@+ probeMuons@-')
)

process.tagMuonsMCMatch = cms.EDProducer("MCTruthDeltaRMatcherNew",
    src = cms.InputTag("tagMuons"),
    matched = cms.InputTag("genParticles"),
    pdgId = cms.vint32(13),
    distMin = cms.double(0.3),
)
process.probeMuonsMCMatch = process.tagMuonsMCMatch.clone(src = "probeMuons")

from MuonAnalysis.TagAndProbe.muon.tag_probe_muon_extraIso_cff import ExtraIsolationVariables

process.load("MuonAnalysis.TagAndProbe.mvaIsoVariables_cff")
from MuonAnalysis.TagAndProbe.mvaIsoVariables_cff import MVAIsoVariablesPlain, MVAIsoVariablesMix
#process.load("MuonAnalysis.TagAndProbe.radialIso_cfi")

process.tpTree = cms.EDAnalyzer("TagProbeFitTreeProducer",
    # choice of tag and probe pairs, and arbitration
    tagProbePairs = cms.InputTag("tpPairs"),
    arbitration   = cms.string("OneProbe"),
    # probe variables: all useful ones
    variables = cms.PSet(
        AllVariables,
        ExtraIsolationVariables,
        MVAIsoVariablesPlain, MVAIsoVariablesMix,
        isoTrk03Abs = cms.InputTag("probeMuonsIsoValueMaps","probeMuonsIsoFromDepsTk"),
        isoTrk03Rel = cms.InputTag("probeMuonsIsoValueMaps","probeMuonsRelIsoFromDepsTk"),
        dxyBS = cms.InputTag("muonDxyPVdzmin","dxyBS"),
        dxyPVdzmin = cms.InputTag("muonDxyPVdzmin","dxyPVdzmin"),
        dzPV = cms.InputTag("muonDxyPVdzmin","dzPV"),
        IP = cms.InputTag("muonSIP","IP"),
        IPError = cms.InputTag("muonSIP","IPError"),
        SIP = cms.InputTag("muonSIP","SIP"),
        #radialIso = cms.InputTag("radialIso"), 
    ),
    flags = cms.PSet(
       TrackQualityFlags,
       MuonIDFlags,
       HighPtTriggerFlags,
    ),
    tagVariables = cms.PSet(
        nVertices   = cms.InputTag("nverticesModule"),
        combRelIso = cms.string("(isolationR03.emEt + isolationR03.hadEt + isolationR03.sumPt)/pt"),
    ),
    tagFlags = cms.PSet(HighPtTriggerFlags),
    pairVariables = cms.PSet(
        nJets30 = cms.InputTag("njets30Module"),
        dz      = cms.string("daughter(0).vz - daughter(1).vz"),
        pt      = cms.string("pt"), # let's do some bump hunt in the T&P too
    ),
    pairFlags = cms.PSet(),
    isMC           = cms.bool(True),
    tagMatches       = cms.InputTag("tagMuonsMCMatch"),
    probeMatches     = cms.InputTag("probeMuonsMCMatch"),
    motherPdgId      = cms.vint32(22, 23),
    makeMCUnbiasTree       = cms.bool(False), 
    checkMotherInUnbiasEff = cms.bool(True),
    allProbes              = cms.InputTag("probeMuons"),
)


process.load('RecoJets.Configuration.RecoPFJets_cff')
##-------------------- Turn-on the FastJet density calculation -----------------------
process.kt6PFJets.doRhoFastjet = True
##-------------------- FastJet density calculation for Iso -----------------------
process.kt6PFJetsForIso = process.kt6PFJets.clone( Rho_EtaMax = cms.double(2.5), Ghost_EtaMax = cms.double(2.5) )


process.load("MuonAnalysis.TagAndProbe.muon.tag_probe_muon_extraIso_cfi")

process.extraProbeVariablesSeq = cms.Sequence(
    process.probeMuonsIsoSequence +
    process.kt6PFJetsForIso * process.computeCorrectedIso + 
    process.mvaIsoVariablesSeq + #* process.radialIso +
    process.muonDxyPVdzmin +
    process.muonSIP 
)
process.tnpSimpleSequence = cms.Sequence(
    process.tagMuons   * process.tagMuonsMCMatch   +
    process.probeMuons * process.probeMuonsMCMatch +
    process.tpPairs    +
    process.nverticesModule +
    process.njets30Module +
    process.extraProbeVariablesSeq +
    process.tpTree
)

#process.fastFilter = cms.EDFilter("MuonSelector", src = cms.InputTag("muons"), cut = cms.string("pt > 15"), filter = cms.bool(True))
process.tagAndProbe = cms.Path( 
    #process.fastFilter                  *
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
process.probeMuonsMCMatchSta = process.tagMuonsMCMatch.clone(src = "probeMuonsSta")

process.tpPairsSta = process.tpPairs.clone(decay = "tagMuons@+ probeMuonsSta@-", cut = '40 < mass < 150')

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
    allProbes     = "probeMuonsSta",
    probeMatches  = "probeMuonsMCMatchSta",
)
process.tpTreeSta.pairVariables.nJets30 = "njets30ModuleSta"
process.njets30ModuleSta = process.njets30Module.clone(pairs = "tpPairsSta")

process.tnpSimpleSequenceSta = cms.Sequence(
    process.tagMuons   * process.tagMuonsMCMatch   +
    process.probeMuonsSta * process.probeMuonsMCMatchSta +
    process.tpPairsSta      +
    process.nverticesModule +
    process.staToTkMatchSequenceZ +
    process.njets30ModuleSta +
    process.tpTreeSta
)

process.tagAndProbeSta = cms.Path( 
    process.muonsSta                       +
    process.patMuonsWithTriggerSequenceSta +
    process.tnpSimpleSequenceSta
)

##    _____     _          ____       _            
##   |  ___|_ _| | _____  |  _ \ __ _| |_ ___  ___ 
##   | |_ / _` | |/ / _ \ | |_) / _` | __/ _ \/ __|
##   |  _| (_| |   <  __/ |  _ < (_| | ||  __/\__ \
##   |_|  \__,_|_|\_\___| |_| \_\__,_|\__\___||___/
##                                                 
##   
process.load("MuonAnalysis.TagAndProbe.fakerate_all_cff")

process.fakeRateJetPlusProbeTree = process.tpTree.clone(
    tagProbePairs = 'jetPlusProbe',
    arbitration   = 'None', 
    tagVariables = process.JetPlusProbeTagVariables,
    tagFlags = cms.PSet(),
    pairVariables = cms.PSet(deltaPhi = cms.string("deltaPhi(daughter(0).phi, daughter(1).phi)")), 
    pairFlags     = cms.PSet(), 
    isMC = False, # MC matches not in place for FR yet
)
process.fakeRateWPlusProbeTree = process.tpTree.clone(
    tagProbePairs = 'wPlusProbe',
    arbitration   = 'None', 
    tagVariables = process.WPlusProbeTagVariables,
    tagFlags = cms.PSet(),
    pairVariables = cms.PSet(), 
    pairFlags     = cms.PSet(SameSign = cms.string('daughter(0).daughter(0).charge == daughter(1).charge')), 
    isMC = False, # MC matches not in place for FR yet
)
process.fakeRateZPlusProbeTree = process.tpTree.clone(
    tagProbePairs = 'zPlusProbe',
    arbitration   = 'None', 
    tagVariables  = process.ZPlusProbeTagVariables,
    tagFlags      = cms.PSet(),
    pairVariables = cms.PSet(), 
    pairFlags     = cms.PSet(), 
    isMC = False, # MC matches not in place for FR yet
)

process.fakeRateJetPlusProbe = cms.Path(
    process.mergedMuons * process.patMuonsWithTriggerSequence +
    process.tagMuons + process.probeMuons + process.extraProbeVariablesSeq + 
    process.jetPlusProbeSequence +
    process.fakeRateJetPlusProbeTree
)
process.fakeRateWPlusProbe = cms.Path(
    process.mergedMuons * process.patMuonsWithTriggerSequence +
    process.tagMuons + process.probeMuons + process.extraProbeVariablesSeq + 
    process.wPlusProbeSequence +
    process.fakeRateWPlusProbeTree
)
process.fakeRateZPlusProbe = cms.Path(
    process.mergedMuons * process.patMuonsWithTriggerSequence +
    process.tagMuons + process.probeMuons + process.extraProbeVariablesSeq + 
    process.zPlusProbeSequence +
    process.fakeRateZPlusProbeTree
)



process.schedule = cms.Schedule(
   process.tagAndProbe, 
   process.tagAndProbeSta, 
   process.fakeRateJetPlusProbe,
   process.fakeRateWPlusProbe,
   process.fakeRateZPlusProbe,
)

process.TFileService = cms.Service("TFileService", fileName = cms.string("tnpZ_MC.root"))

import FWCore.ParameterSet.Config as cms

import subprocess
#dataSummary = open('dataSummary.txt', 'w')

process = cms.Process("TagProbe")

process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.source = cms.Source("PoolSource", 
    fileNames = cms.untracked.vstring(),
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )    

process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")

import os
if "CMSSW_7_4_" in os.environ['CMSSW_VERSION']:

    #run 251168
    process.GlobalTag.globaltag = cms.string('GR_P_V56::All')
    sourcefilesfolder = "/store/data/Run2015B/SingleMuon/AOD/PromptReco-v1/000/251/168/00000"
    files = subprocess.check_output([ "/afs/cern.ch/project/eos/installation/0.3.15/bin/eos.select", "ls", sourcefilesfolder ])
    process.source.fileNames = [ sourcefilesfolder+"/"+f for f in files.split() ]

    #run 251244
    sourcefilesfolder = "/store/data/Run2015B/SingleMuon/AOD/PromptReco-v1/000/251/244/00000"
    files = subprocess.check_output([ "/afs/cern.ch/project/eos/installation/0.3.15/bin/eos.select", "ls", sourcefilesfolder ])
    process.source.fileNames.extend( [ sourcefilesfolder+"/"+f for f in files.split() ] )

    #run 251251
    sourcefilesfolder = "/store/data/Run2015B/SingleMuon/AOD/PromptReco-v1/000/251/251/00000"
    files = subprocess.check_output([ "/afs/cern.ch/project/eos/installation/0.3.15/bin/eos.select", "ls", sourcefilesfolder ])
    process.source.fileNames.extend( [ sourcefilesfolder+"/"+f for f in files.split() ] )

    #run 251252
    sourcefilesfolder = "/store/data/Run2015B/SingleMuon/AOD/PromptReco-v1/000/251/252/00000"
    files = subprocess.check_output([ "/afs/cern.ch/project/eos/installation/0.3.15/bin/eos.select", "ls", sourcefilesfolder ])
    process.source.fileNames.extend( [ sourcefilesfolder+"/"+f for f in files.split() ] )

    # to add following runs: 251491, 251493, 251496, ..., 251500 
    print process.source.fileNames
    #print process.source.fileNames, dataSummary

else: raise RuntimeError, "Unknown CMSSW version %s" % os.environ['CMSSW_VERSION']

## SELECT WHAT DATASET YOU'RE RUNNING ON
TRIGGER="SingleMu"
#TRIGGER="DoubleMu"

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

process.load("HLTrigger.HLTfilters.triggerResultsFilter_cfi")


if TRIGGER == "SingleMu":
    process.triggerResultsFilter.triggerConditions = cms.vstring( 'HLT_IsoMu24_eta2p1_v*', 'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*' )
elif TRIGGER == "DoubleMu":
    process.triggerResultsFilter.triggerConditions = cms.vstring( 'HLT_Mu8_v*', 'HLT_Mu17_v*', 'HLT_Mu17_TkMu8_NoDZ_v*', 'HLT_Mu13_Mu8_NoDZ_v*' )
else:
    raise RuntimeError, "TRIGGER must be 'SingleMu' or 'DoubleMu'"

process.triggerResultsFilter.l1tResults = "gtDigis"
process.triggerResultsFilter.throw = False
process.triggerResultsFilter.hltResults = cms.InputTag("TriggerResults","","HLT")

#decomment when you have it
#process.triggerResultsFilterFake = process.triggerResultsFilter.clone(
#    triggerConditions = cms.vstring( 'HLT_Mu40_v*', 'HLT_Mu5_v*', 'HLT_Mu12_v*', 'HLT_Mu24_v*')
#)

process.fastFilter     = cms.Sequence(process.goodVertexFilter + process.noScraping + process.triggerResultsFilter)
#process.fastFilterFake = cms.Sequence(process.goodVertexFilter + process.noScraping + process.triggerResultsFilterFake)

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
    muonsCut     = cms.string("pt > 3 && track.isNonnull"),
    caloMuonsCut = cms.string("pt > 3"),
    tracksCut    = cms.string("pt > 3"),
)

## ==== Trigger matching
process.load("MuonAnalysis.MuonAssociators.patMuonsWithTrigger_cff")
## with some customization
process.muonMatchHLTL2.maxDeltaR = 0.3 # Zoltan tuning - it was 0.5
process.muonMatchHLTL3.maxDeltaR = 0.1
from MuonAnalysis.MuonAssociators.patMuonsWithTrigger_cff import *
changeRecoMuonInput(process, "mergedMuons")
#useExtendedL1Match(process)
#addHLTL1Passthrough(process)


from MuonAnalysis.TagAndProbe.common_variables_cff import *
process.load("MuonAnalysis.TagAndProbe.common_modules_cff")

process.tagMuons = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("patMuonsWithTrigger"),
    cut = cms.string("pt > 15 && "+MuonIDFlags.Tight2012.value()+
                     " && !triggerObjectMatchesByCollection('hltL3MuonCandidates').empty()"+
                     " && pfIsolationR04().sumChargedHadronPt/pt < 0.2"),
)
if TRIGGER == "DoubleMu":
    process.tagMuons.cut = ("pt > 6 && (isGlobalMuon || isTrackerMuon) && isPFMuon "+
                            " && !triggerObjectMatchesByCollection('hltL3MuonCandidates').empty()"+
                            " && pfIsolationR04().sumChargedHadronPt/pt < 0.2")

process.oneTag  = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("tagMuons"), minNumber = cms.uint32(1))

process.probeMuons = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("patMuonsWithTrigger"),
    cut = cms.string("track.isNonnull"),  # no real cut now
)

process.tpPairs = cms.EDProducer("CandViewShallowCloneCombiner",
    #cut = cms.string('60 < mass < 140 && abs(daughter(0).vz - daughter(1).vz) < 4'),
    cut = cms.string('60 < mass && abs(daughter(0).vz - daughter(1).vz) < 4'),
    decay = cms.string('tagMuons@+ probeMuons@-')
)
process.onePair = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("tpPairs"), minNumber = cms.uint32(1))

from MuonAnalysis.TagAndProbe.muon.tag_probe_muon_extraIso_cff import ExtraIsolationVariables

process.load("MuonAnalysis.TagAndProbe.mvaIsoVariables_cff")
from MuonAnalysis.TagAndProbe.mvaIsoVariables_cff import MVAIsoVariablesPlain, MVAIsoVariablesPlainTag
process.load("MuonAnalysis.TagAndProbe.radialIso_cfi")

process.tpTree = cms.EDAnalyzer("TagProbeFitTreeProducer",
    # choice of tag and probe pairs, and arbitration
    tagProbePairs = cms.InputTag("tpPairs"),
    arbitration   = cms.string("None"),
    # probe variables: all useful ones
    variables = cms.PSet(
        AllVariables,
        ExtraIsolationVariables,
        MVAIsoVariablesPlain, 
        isoTrk03Abs = cms.InputTag("probeMuonsIsoValueMaps","probeMuonsIsoFromDepsTk"),
        isoTrk03Rel = cms.InputTag("probeMuonsIsoValueMaps","probeMuonsRelIsoFromDepsTk"),
        dxyBS = cms.InputTag("muonDxyPVdzmin","dxyBS"),
        dxyPVdzmin = cms.InputTag("muonDxyPVdzmin","dxyPVdzmin"),
        dzPV = cms.InputTag("muonDxyPVdzmin","dzPV"),
        radialIso = cms.InputTag("radialIso"), 
        nSplitTk  = cms.InputTag("splitTrackTagger"),
    ),
    flags = cms.PSet(
       TrackQualityFlags,
       MuonIDFlags,
       HighPtTriggerFlags,
       HighPtTriggerFlagsDebug,
    ),
    tagVariables = cms.PSet(
     #   TriggerVariables, 
     #   MVAIsoVariablesPlainTag, 
     #   pt = cms.string("pt"),
     #   eta = cms.string("eta"),
     #   phi = cms.string("phi"),
     #   nVertices   = cms.InputTag("nverticesModule"),
     #   combRelIso = cms.string("(isolationR03.emEt + isolationR03.hadEt + isolationR03.sumPt)/pt"),
     #   chargedHadIso04 = cms.string("pfIsolationR04().sumChargedHadronPt"),
     #   neutralHadIso04 = cms.string("pfIsolationR04().sumNeutralHadronEt"),
     #   photonIso04 = cms.string("pfIsolationR04().sumPhotonEt"),
     #   combRelIsoPF04dBeta = IsolationVariables.combRelIsoPF04dBeta,
     #   combRelIsoPF03dBeta = IsolationVariables.combRelIsoPF03dBeta,
     #   dzPV = cms.InputTag("muonDxyPVdzminTags","dzPV"),
        AllVariables,
        ExtraIsolationVariables,
        MVAIsoVariablesPlain, 
        isoTrk03Abs = cms.InputTag("probeMuonsIsoValueMaps","probeMuonsIsoFromDepsTk"),
        isoTrk03Rel = cms.InputTag("probeMuonsIsoValueMaps","probeMuonsRelIsoFromDepsTk"),
        dxyBS = cms.InputTag("muonDxyPVdzminTags","dxyBS"),
        dxyPVdzmin = cms.InputTag("muonDxyPVdzminTags","dxyPVdzmin"),
        dzPV = cms.InputTag("muonDxyPVdzminTags","dzPV"),
        radialIso = cms.InputTag("radialIso"), 
        nSplitTk  = cms.InputTag("splitTrackTagger"),
        #l1rate = cms.InputTag("l1rate"),
        #bx     = cms.InputTag("l1rate","bx"),
        #mu17ps = cms.InputTag("l1hltprescale","HLTMu17TotalPrescale"), 
        #mu8ps  = cms.InputTag("l1hltprescale","HLTMu8TotalPrescale"), 
    ),
    tagFlags = cms.PSet(HighPtTriggerFlags,HighPtTriggerFlagsDebug),
    pairVariables = cms.PSet(
        nJets30 = cms.InputTag("njets30Module"),
        dz      = cms.string("daughter(0).vz - daughter(1).vz"),
        pt      = cms.string("pt"), 
        rapidity = cms.string("rapidity"),
        deltaR   = cms.string("deltaR(daughter(0).eta, daughter(0).phi, daughter(1).eta, daughter(1).phi)"), 
        probeMultiplicity = cms.InputTag("probeMultiplicity"),
        probeMultiplicity_TMGM = cms.InputTag("probeMultiplicityTMGM"),
        probeMultiplicity_Pt10_M60140 = cms.InputTag("probeMultiplicityPt10M60140"),
        ## New TuneP variables
        newTuneP_probe_pt            = cms.InputTag("newTunePVals", "pt"),
        newTuneP_probe_sigmaPtOverPt = cms.InputTag("newTunePVals", "ptRelError"),
        newTuneP_probe_trackType     = cms.InputTag("newTunePVals", "trackType"),
        newTuneP_mass                = cms.InputTag("newTunePVals", "mass"),
    ),
    pairFlags = cms.PSet(
        BestZ = cms.InputTag("bestPairByZMass"),
    ),
    isMC           = cms.bool(False),
    addRunLumiInfo = cms.bool(True),
)
if TRIGGER == "DoubleMu":
    for K,F in MuonIDFlags.parameters_().iteritems():
        setattr(process.tpTree.tagFlags, K, F)


process.load("MuonAnalysis.TagAndProbe.muon.tag_probe_muon_extraIso_cfi")

process.extraProbeVariablesSeq = cms.Sequence(
    process.probeMuonsIsoSequence +
    process.computeCorrectedIso + 
    process.mvaIsoVariablesSeq * process.mvaIsoVariablesTag * process.radialIso +
    process.splitTrackTagger +
    process.muonDxyPVdzmin 
)

process.tnpSimpleSequence = cms.Sequence(
    process.tagMuons +
    process.oneTag     +
    process.probeMuons +
    process.tpPairs    +
    process.onePair    +
    process.nverticesModule +
    process.njets30Module +
    process.extraProbeVariablesSeq +
    process.probeMultiplicities + 
    #process.l1rate +
    #process.l1hltprescale + 
    process.bestPairByZMass + 
    process.newTunePVals +
    process.muonDxyPVdzminTags +
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

process.tpPairsSta = process.tpPairs.clone(decay = "tagMuons@+ probeMuonsSta@-", cut = '40 < mass < 150')

process.onePairSta = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("tpPairsSta"), minNumber = cms.uint32(1))

process.staToTkMatch.maxDeltaR     = 0.3
process.staToTkMatch.maxDeltaPtRel = 2.
process.staToTkMatchNoZ.maxDeltaR     = 0.3
process.staToTkMatchNoZ.maxDeltaPtRel = 2.

process.load("MuonAnalysis.TagAndProbe.tracking_reco_info_cff")

process.tpTreeSta = process.tpTree.clone(
    tagProbePairs = "tpPairsSta",
    arbitration   = "OneProbe",
    variables = cms.PSet(
        KinematicVariables, 
        StaOnlyVariables,
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
        Tk  = cms.string("track.isNonnull"),
        StaTkSameCharge = cms.string("outerTrack.isNonnull && innerTrack.isNonnull && (outerTrack.charge == innerTrack.charge)"),
    ),
    tagVariables = cms.PSet(
        pt = cms.string("pt"),
        eta = cms.string("eta"),
        phi = cms.string("phi"),
        nVertices = cms.InputTag("nverticesModule"),
        combRelIso = cms.string("(isolationR03.emEt + isolationR03.hadEt + isolationR03.sumPt)/pt"),
        chargedHadIso04 = cms.string("pfIsolationR04().sumChargedHadronPt"),
        neutralHadIso04 = cms.string("pfIsolationR04().sumNeutralHadronEt"),
        photonIso04 = cms.string("pfIsolationR04().sumPhotonEt"),
        combRelIsoPF04dBeta = IsolationVariables.combRelIsoPF04dBeta,
        #l1rate = cms.InputTag("l1rate"),
        #bx     = cms.InputTag("l1rate","bx"),
    ),
    pairVariables = cms.PSet(
        nJets30 = cms.InputTag("njets30ModuleSta"),
        dz      = cms.string("daughter(0).vz - daughter(1).vz"),
        pt      = cms.string("pt"), 
        rapidity = cms.string("rapidity"),
        deltaR   = cms.string("deltaR(daughter(0).eta, daughter(0).phi, daughter(1).eta, daughter(1).phi)"), 
    ),
    pairFlags = cms.PSet(),
)
process.njets30ModuleSta = process.njets30Module.clone(pairs = "tpPairsSta")

process.tnpSimpleSequenceSta = cms.Sequence(
    process.tagMuons +
    process.oneTag     +
    process.probeMuonsSta   +
    process.tpPairsSta      +
    process.onePairSta      +
    process.nverticesModule +
    process.staToTkMatchSequenceZ +
    process.njets30ModuleSta +
    #process.l1rate +
    process.tpTreeSta
)

## Add extra RECO-level info
if False:
    process.tnpSimpleSequenceSta.replace(process.tpTreeSta, process.tkClusterInfo+process.tpTreeSta)
    process.tpTreeSta.tagVariables.nClustersStrip = cms.InputTag("tkClusterInfo","siStripClusterCount")
    process.tpTreeSta.tagVariables.nClustersPixel = cms.InputTag("tkClusterInfo","siPixelClusterCount")
    process.tnpSimpleSequenceSta.replace(process.tpTreeSta, process.tkLogErrors+process.tpTreeSta)
    process.tpTreeSta.tagVariables.nLogErrFirst = cms.InputTag("tkLogErrors","firstStep")
    process.tpTreeSta.tagVariables.nLogErrPix   = cms.InputTag("tkLogErrors","pixelSteps")
    process.tpTreeSta.tagVariables.nLogErrAny   = cms.InputTag("tkLogErrors","anyStep")

if True: 
    process.tracksNoMuonSeeded = cms.EDFilter("TrackSelector",
      src = cms.InputTag("generalTracks"),
      cut = cms.string(" || ".join("isAlgoInMask('%s')" % a for a in [
                    'initialStep', 'lowPtTripletStep', 'pixelPairStep', 'detachedTripletStep',
                    'mixedTripletStep', 'pixelLessStep', 'tobTecStep', 'jetCoreRegionalStep' ] ) )
    )
    process.pCutTracks0 = process.pCutTracks.clone(src = 'tracksNoMuonSeeded')
    process.tkTracks0 = process.tkTracks.clone(src = 'pCutTracks0')
    process.tkTracksNoZ0 = process.tkTracksNoZ.clone(src = 'tkTracks0')
    process.preTkMatchSequenceZ.replace(
            process.tkTracksNoZ, process.tkTracksNoZ +
            process.tracksNoMuonSeeded + process.pCutTracks0 + process.tkTracks0 + process.tkTracksNoZ0)
    process.staToTkMatch0 = process.staToTkMatch.clone(matched = 'tkTracks0')
    process.staToTkMatchNoZ0 = process.staToTkMatchNoZ.clone(matched = 'tkTracksNoZ0')
    process.staToTkMatchSequenceZ.replace( process.staToTkMatch, process.staToTkMatch + process.staToTkMatch0 )
    process.staToTkMatchSequenceZ.replace( process.staToTkMatchNoZ, process.staToTkMatchNoZ + process.staToTkMatchNoZ0 )
    process.tpTreeSta.variables.tk0_deltaR     = cms.InputTag("staToTkMatch0","deltaR")
    process.tpTreeSta.variables.tk0_deltaEta   = cms.InputTag("staToTkMatch0","deltaEta")
    process.tpTreeSta.variables.tk0_deltaR_NoZ   = cms.InputTag("staToTkMatchNoZ0","deltaR")
    process.tpTreeSta.variables.tk0_deltaEta_NoZ = cms.InputTag("staToTkMatchNoZ0","deltaEta")

process.tagAndProbeSta = cms.Path( 
    process.fastFilter +
    process.muonsSta                       +
    process.patMuonsWithTriggerSequenceSta +
    process.tnpSimpleSequenceSta
)


if True: # turn on for tracking efficiency using L1 seeds
    process.probeL1 = cms.EDFilter("CandViewSelector",
        src = cms.InputTag("l1extraParticles"),
        cut = cms.string("pt >= 5 && abs(eta) < 2.4"),
    )
    process.tpPairsTkL1 = process.tpPairs.clone(decay = "tagMuons@+ probeL1@-", cut = 'mass > 30')
    process.l1ToTkMatch    = process.staToTkMatch.clone(src = "probeL1", srcTrack="none")
    process.l1ToTkMatchNoZ = process.staToTkMatchNoZ.clone(src = "probeL1", srcTrack="none")
    process.l1ToTkMatch0    = process.staToTkMatch0.clone(src = "probeL1", srcTrack="none")
    process.l1ToTkMatchNoZ0 = process.staToTkMatchNoZ0.clone(src = "probeL1", srcTrack="none")
    process.tpTreeL1 = process.tpTreeSta.clone(
        tagProbePairs = "tpPairsTkL1",
        arbitration   = "OneProbe",
        variables = cms.PSet(
            KinematicVariables, 
            bx = cms.string("bx"),
            quality = cms.string("gmtMuonCand.quality"),
            ## track matching variables
            tk_deltaR     = cms.InputTag("l1ToTkMatch","deltaR"),
            tk_deltaEta   = cms.InputTag("l1ToTkMatch","deltaEta"),
            tk_deltaR_NoZ   = cms.InputTag("l1ToTkMatchNoZ","deltaR"),
            tk_deltaEta_NoZ = cms.InputTag("l1ToTkMatchNoZ","deltaEta"),
            ## track matching variables (early general tracks)
            tk0_deltaR     = cms.InputTag("l1ToTkMatch0","deltaR"),
            tk0_deltaEta   = cms.InputTag("l1ToTkMatch0","deltaEta"),
            tk0_deltaR_NoZ   = cms.InputTag("l1ToTkMatchNoZ0","deltaR"),
            tk0_deltaEta_NoZ = cms.InputTag("l1ToTkMatchNoZ0","deltaEta"),
        ),
        flags = cms.PSet(
        ),
        tagVariables = cms.PSet(
            pt = cms.string("pt"),
            eta = cms.string("eta"),
            phi = cms.string("phi"),
            nVertices   = cms.InputTag("nverticesModule"),
            combRelIso = cms.string("(isolationR03.emEt + isolationR03.hadEt + isolationR03.sumPt)/pt"),
            chargedHadIso04 = cms.string("pfIsolationR04().sumChargedHadronPt"),
            neutralHadIso04 = cms.string("pfIsolationR04().sumNeutralHadronEt"),
            photonIso04 = cms.string("pfIsolationR04().sumPhotonEt"),
            combRelIsoPF04dBeta = IsolationVariables.combRelIsoPF04dBeta,
        ),
        pairVariables = cms.PSet(
            #nJets30 = cms.InputTag("njets30ModuleSta"),
            pt      = cms.string("pt"),
            rapidity = cms.string("rapidity"),
            deltaR   = cms.string("deltaR(daughter(0).eta, daughter(0).phi, daughter(1).eta, daughter(1).phi)"), 
        ),
        pairFlags = cms.PSet(),
        allProbes     = cms.InputTag("probeL1"),
    )
    process.tagAndProbeTkL1 = cms.Path(
        process.fastFilter +
        process.probeL1 +
        process.tpPairsTkL1 +
        process.preTkMatchSequenceZ +
        process.l1ToTkMatch + process.l1ToTkMatchNoZ +
        process.l1ToTkMatch0 + process.l1ToTkMatchNoZ0 +
        process.tpTreeL1
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
)
process.fakeRateWPlusProbeTree = process.tpTree.clone(
    tagProbePairs = 'wPlusProbe',
    arbitration   = 'None', 
    tagVariables = process.WPlusProbeTagVariables,
    tagFlags = cms.PSet(),
    pairVariables = cms.PSet(), 
    pairFlags     = cms.PSet(SameSign = cms.string('daughter(0).daughter(0).charge == daughter(1).charge')), 
)
process.fakeRateZPlusProbeTree = process.tpTree.clone(
    tagProbePairs = 'zPlusProbe',
    arbitration   = 'None', 
    tagVariables  = process.ZPlusProbeTagVariables,
    tagFlags      = cms.PSet(),
    pairVariables = cms.PSet(), 
    pairFlags     = cms.PSet(), 
)

process.fakeRateJetPlusProbe = cms.Path(
    #process.fastFilterFake +
    process.mergedMuons * process.patMuonsWithTriggerSequence +
    process.tagMuons + process.probeMuons + process.extraProbeVariablesSeq + 
    process.jetPlusProbeSequence +
    process.fakeRateJetPlusProbeTree
)
process.fakeRateWPlusProbe = cms.Path(
    process.fastFilter +
    process.mergedMuons * process.patMuonsWithTriggerSequence +
    process.tagMuons + process.probeMuons + process.extraProbeVariablesSeq + 
    process.wPlusProbeSequence +
    process.fakeRateWPlusProbeTree
)
process.fakeRateZPlusProbe = cms.Path(
    process.fastFilter +
    process.mergedMuons * process.patMuonsWithTriggerSequence +
    process.tagMuons + process.probeMuons + process.extraProbeVariablesSeq + 
    process.zPlusProbeSequence +
    process.fakeRateZPlusProbeTree
)

process.schedule = cms.Schedule(
   process.tagAndProbe, 
   process.tagAndProbeSta, 
   process.tagAndProbeTkL1
)

process.RandomNumberGeneratorService.tkTracksNoZ = cms.PSet( initialSeed = cms.untracked.uint32(81) )
process.RandomNumberGeneratorService.tkTracksNoZ0 = cms.PSet( initialSeed = cms.untracked.uint32(81) )


if TRIGGER == "SingleMu": 
    process.schedule.extend([
       process.fakeRateJetPlusProbe,
       process.fakeRateWPlusProbe,
       process.fakeRateZPlusProbe,
    ])

process.TFileService = cms.Service("TFileService", fileName = cms.string("tnpZ_Data_PromptReco.root"))

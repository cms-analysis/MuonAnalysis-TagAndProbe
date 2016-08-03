import FWCore.ParameterSet.Config as cms

import subprocess

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

# process.Tracer = cms.Service('Tracer')
import os
if   "CMSSW_5_3_" in os.environ['CMSSW_VERSION']:
    process.GlobalTag.globaltag = cms.string('GR_R_53_V7::All')
    process.source.fileNames = [
        '/store/data/Run2012C/MuOnia/AOD/PromptReco-v1/000/198/208/B0C8994D-2CC7-E111-B04E-0025901D6288.root',
        '/store/data/Run2012C/MuOnia/AOD/PromptReco-v1/000/198/208/76FB8276-3AC7-E111-A52E-001D09F297EF.root',
        '/store/data/Run2012C/MuOnia/AOD/PromptReco-v1/000/198/208/76E299E2-2FC7-E111-B2CC-001D09F28F25.root',
        '/store/data/Run2012C/MuOnia/AOD/PromptReco-v1/000/198/208/5C0AD4DE-4AC7-E111-9C14-001D09F23A20.root',
        '/store/data/Run2012C/MuOnia/AOD/PromptReco-v1/000/198/208/4E4682A9-30C7-E111-9A8C-003048F1C420.root',
        '/store/data/Run2012C/MuOnia/AOD/PromptReco-v1/000/198/208/3669603F-2FC7-E111-8935-003048D2BE06.root',
        '/store/data/Run2012C/MuOnia/AOD/PromptReco-v1/000/198/208/0E9130E2-2FC7-E111-9494-001D09F26509.root',
        '/store/data/Run2012C/MuOnia/AOD/PromptReco-v1/000/198/208/02CA3FBC-2DC7-E111-9CA1-003048D2BB90.root',
    ]
elif "CMSSW_5_2_" in os.environ['CMSSW_VERSION']:
    process.GlobalTag.globaltag = cms.string('GR_P_V39_AN1::All')
    process.source.fileNames = [
        '/store/data/Run2012A/MuOnia/AOD/PromptReco-v1/000/191/226/645E9BE9-FC87-E111-9D30-5404A63886C5.root',
        '/store/data/Run2012A/MuOnia/AOD/PromptReco-v1/000/191/226/60C36C17-0188-E111-8B3D-5404A638868F.root',
        '/store/data/Run2012A/MuOnia/AOD/PromptReco-v1/000/191/226/5E384785-F087-E111-A8CB-BCAEC518FF80.root',
        '/store/data/Run2012A/MuOnia/AOD/PromptReco-v1/000/191/226/5C605E9F-F887-E111-A540-00237DDC5CB0.root',
        '/store/data/Run2012A/MuOnia/AOD/PromptReco-v1/000/191/226/5C0DCA92-CE87-E111-A539-5404A63886C6.root',
        '/store/data/Run2012A/MuOnia/AOD/PromptReco-v1/000/191/226/5A856F46-1888-E111-A3A5-BCAEC5329700.root',
        '/store/data/Run2012A/MuOnia/AOD/PromptReco-v1/000/191/226/5271187F-DF87-E111-ADDB-BCAEC518FF7C.root',
    ]
elif "CMSSW_7_4_" in os.environ['CMSSW_VERSION']:
    process.GlobalTag.globaltag = cms.string('GR_P_V56::All')
    process.source.fileNames = [
        '/store/data/Run2015B/Charmonium/AOD/PromptReco-v1/000/251/164/00000/56BA3AC3-A226-E511-98BB-02163E0133B5.root',
	'/store/data/Run2015B/Charmonium/AOD/PromptReco-v1/000/251/167/00000/FED37D13-A826-E511-8641-02163E01386E.root',
	'/store/data/Run2015B/Charmonium/AOD/PromptReco-v1/000/251/168/00000/1602316B-CF26-E511-BCBA-02163E011BF3.root',
	'/store/data/Run2015B/Charmonium/AOD/PromptReco-v1/000/251/168/00000/86456DCF-D026-E511-A8E2-02163E011D37.root',
        '/store/data/Run2015B/Charmonium/AOD/PromptReco-v1/000/251/162/00000/E855BF93-4227-E511-8207-02163E011976.root',
        '/store/data/Run2015B/Charmonium/AOD/PromptReco-v1/000/251/244/00000/4ADDEB99-6727-E511-BEF2-02163E011955.root',
        '/store/data/Run2015B/Charmonium/AOD/PromptReco-v1/000/251/244/00000/C6DB3E01-8327-E511-B696-02163E0139CF.root',
        '/store/data/Run2015B/Charmonium/AOD/PromptReco-v1/000/251/244/00000/FC4B5FE2-8A27-E511-A0DD-02163E014729.root',
        '/store/data/Run2015B/Charmonium/AOD/PromptReco-v1/000/251/251/00000/7097B2BC-8E27-E511-80E1-02163E0138B3.root',
        '/store/data/Run2015B/Charmonium/AOD/PromptReco-v1/000/251/252/00000/0067F5A4-9A27-E511-B904-02163E01267F.root',
        '/store/data/Run2015B/Charmonium/AOD/PromptReco-v1/000/251/252/00000/642ACBD2-A127-E511-A59E-02163E0123F1.root',
        '/store/data/Run2015B/Charmonium/AOD/PromptReco-v1/000/251/252/00000/F24C3492-9627-E511-AEC5-02163E011C7F.root'
    ]
elif "CMSSW_7_6_" in os.environ['CMSSW_VERSION']:
    process.GlobalTag.globaltag = cms.string('76X_dataRun2_v15')
    process.source.fileNames = [
        '/store/data/Run2015D/Charmonium/AOD/16Dec2015-v1/50000/0013A0EA-94AE-E511-9B84-001E67398683.root',
        '/store/data/Run2015D/Charmonium/AOD/16Dec2015-v1/50000/02A9F67F-A2AE-E511-939A-0CC47A4C8EB6.root',
        '/store/data/Run2015D/Charmonium/AOD/16Dec2015-v1/50000/02D39573-BDAE-E511-AF22-00266CF9AF00.root',
        '/store/data/Run2015D/Charmonium/AOD/16Dec2015-v1/50000/0639B2A3-A7AE-E511-887E-A0000420FE80.root',
        '/store/data/Run2015D/Charmonium/AOD/16Dec2015-v1/50000/06C00544-AAAE-E511-831D-00A0D1EE2BBC.root',
        '/store/data/Run2015D/Charmonium/AOD/16Dec2015-v1/50000/06DCCB30-C7AE-E511-A318-0025905A6092.root',
        '/store/data/Run2015D/Charmonium/AOD/16Dec2015-v1/50000/0ADF8F02-94AE-E511-96F0-549F358EB72E.root',
        '/store/data/Run2015D/Charmonium/AOD/16Dec2015-v1/50000/0C78D007-95AE-E511-87A9-0CC47A4DEDCA.root',
        '/store/data/Run2015D/Charmonium/AOD/16Dec2015-v1/50000/10438105-95AE-E511-ADCD-0090FAA59634.root',
        '/store/data/Run2015D/Charmonium/AOD/16Dec2015-v1/50000/1077A155-AAAE-E511-A1C7-0CC47A4DEF3E.root',
    ]
elif "CMSSW_8_0_" in os.environ['CMSSW_VERSION']:
    process.GlobalTag.globaltag = cms.string('80X_dataRun2_Prompt_v9')
    process.source.fileNames = [
        'root://cms-xrd-global.cern.ch//store/data/Run2016B/Charmonium/AOD/PromptReco-v2/000/274/959/00000/0645C458-9731-E611-BE85-02163E01466B.root',
        'root://cms-xrd-global.cern.ch//store/data/Run2016B/Charmonium/AOD/PromptReco-v2/000/274/959/00000/2AFCDB80-9A31-E611-A6AF-02163E01347E.root',
        'root://cms-xrd-global.cern.ch//store/data/Run2016B/Charmonium/AOD/PromptReco-v2/000/274/959/00000/5427165D-A231-E611-BFF6-02163E01273F.root',
        'root://cms-xrd-global.cern.ch//store/data/Run2016B/Charmonium/AOD/PromptReco-v2/000/274/959/00000/6848311F-AE31-E611-B032-02163E011F8C.root',
        'root://cms-xrd-global.cern.ch//store/data/Run2016B/Charmonium/AOD/PromptReco-v2/000/274/959/00000/7084CDFE-9831-E611-A940-02163E0146A0.root',
        'root://cms-xrd-global.cern.ch//store/data/Run2016B/Charmonium/AOD/PromptReco-v2/000/274/959/00000/740BC96E-9831-E611-9203-02163E011F87.root',
        'root://cms-xrd-global.cern.ch//store/data/Run2016B/Charmonium/AOD/PromptReco-v2/000/274/959/00000/82A366D4-9931-E611-BB5E-02163E011DC5.root',
        'root://cms-xrd-global.cern.ch//store/data/Run2016B/Charmonium/AOD/PromptReco-v2/000/274/959/00000/B8CEAB02-9831-E611-9E77-02163E0134BC.root',
        'root://cms-xrd-global.cern.ch//store/data/Run2016B/Charmonium/AOD/PromptReco-v2/000/274/959/00000/D4805CD1-9531-E611-B434-02163E01339C.root',
        'root://cms-xrd-global.cern.ch//store/data/Run2016B/Charmonium/AOD/PromptReco-v2/000/274/959/00000/D604489F-9C31-E611-A0A1-02163E0135A8.root',
        'root://cms-xrd-global.cern.ch//store/data/Run2016B/Charmonium/AOD/PromptReco-v2/000/274/959/00000/EA753A39-9531-E611-8F4F-02163E014775.root',
        'root://cms-xrd-global.cern.ch//store/data/Run2016B/Charmonium/AOD/PromptReco-v2/000/274/959/00000/EC93E1BC-9F31-E611-8538-02163E0146B6.root',
        'root://cms-xrd-global.cern.ch//store/data/Run2016B/Charmonium/AOD/PromptReco-v2/000/274/959/00000/F2E77417-9431-E611-A181-02163E011EB0.root',
        'root://cms-xrd-global.cern.ch//store/data/Run2016B/Charmonium/AOD/PromptReco-v2/000/274/959/00000/F460DF68-9E31-E611-96DF-02163E011F4D.root',
    ]

else: raise RuntimeError, "Unknown CMSSW version %s" % os.environ['CMSSW_VERSION']



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
process.triggerResultsFilter.triggerConditions = cms.vstring( 'HLT_Mu*_L2Mu*' )
process.triggerResultsFilter.l1tResults = ''
process.triggerResultsFilter.throw = True
process.triggerResultsFilter.hltResults = cms.InputTag( "TriggerResults", "", "HLT" )
process.HLTMu   = process.triggerResultsFilter.clone(triggerConditions = [ 'HLT_Mu*_L2Mu*' ])
process.HLTBoth = process.triggerResultsFilter.clone(triggerConditions = [ 'HLT_Mu*_L2Mu*', 'HLT_Mu*_Track*_Jpsi*' ])

process.fastFilter = cms.Sequence(process.goodVertexFilter + process.noScraping)
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
    muonsCut     = cms.string("pt > 2 && track.isNonnull"),
    caloMuonsCut = cms.string("pt > 2"),
    tracksCut    = cms.string("pt > 2"),
)

## ==== Trigger matching
process.load("MuonAnalysis.MuonAssociators.patMuonsWithTrigger_cff")
## with some customization
from MuonAnalysis.MuonAssociators.patMuonsWithTrigger_cff import *
changeRecoMuonInput(process, "mergedMuons")
#useExtendedL1Match(process)
#addHLTL1Passthrough(process)


from MuonAnalysis.TagAndProbe.common_variables_cff import *
process.load("MuonAnalysis.TagAndProbe.common_modules_cff")

process.tagMuons = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("patMuonsWithTrigger"),
    cut = cms.string("(isGlobalMuon || numberOfMatchedStations > 1) && pt > 5 && !triggerObjectMatchesByCollection('hltL3MuonCandidates').empty()"),
)

process.oneTag  = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("tagMuons"), minNumber = cms.uint32(1))

process.probeMuons = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("patMuonsWithTrigger"),
    cut = cms.string("track.isNonnull && (!triggerObjectMatchesByCollection('hltTracksIter').empty() || !triggerObjectMatchesByCollection('hltMuTrackJpsiEffCtfTrackCands').empty() || !triggerObjectMatchesByCollection('hltMuTrackJpsiCtfTrackCands').empty() || !triggerObjectMatchesByCollection('hltL2MuonCandidates').empty())"),
)

process.tpPairs = cms.EDProducer("CandViewShallowCloneCombiner",
    cut = cms.string('2.8 < mass < 3.4 && abs(daughter(0).vz - daughter(1).vz) < 1'),
    decay = cms.string('tagMuons@+ probeMuons@-')
)
process.onePair = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("tpPairs"), minNumber = cms.uint32(1))

from MuonAnalysis.TagAndProbe.muon.tag_probe_muon_extraIso_cff import ExtraIsolationVariables

process.tpTree = cms.EDAnalyzer("TagProbeFitTreeProducer",
    # choice of tag and probe pairs, and arbitration
    tagProbePairs = cms.InputTag("tpPairs"),
    arbitration   = cms.string("None"),
    # probe variables: all useful ones
    variables = cms.PSet(
        AllVariables,
        ExtraIsolationVariables,
        dxyBS = cms.InputTag("muonDxyPVdzmin","dxyBS"),
        dxyPVdzmin = cms.InputTag("muonDxyPVdzmin","dxyPVdzmin"),
        dzPV = cms.InputTag("muonDxyPVdzmin","dzPV"),
        nSplitTk  = cms.InputTag("splitTrackTagger"),
    ),
    flags = cms.PSet(
       TrackQualityFlags,
       MuonIDFlags,
       HighPtTriggerFlags,
       HighPtTriggerFlagsDebug,
       LowPtTriggerFlagsPhysics,
       LowPtTriggerFlagsEfficienciesProbe,
       Acc_JPsi = cms.string("(abs(eta) <= 1.3 && pt > 3.3) || (1.3 < abs(eta) <= 2.2 && p > 2.9) || (2.2 < abs(eta) <= 2.4  && pt > 0.8)"),
    ),
    tagVariables = cms.PSet(
        pt  = cms.string('pt'),
        eta = cms.string('eta'),
        phi = cms.string('phi'),
        nVertices = cms.InputTag("nverticesModule"),
        l1rate = cms.InputTag("l1rate"),
        bx     = cms.InputTag("l1rate","bx"),
    ),
    tagFlags     = cms.PSet(
        HighPtTriggerFlags,
        HighPtTriggerFlagsDebug,
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
        dz            = cms.string("daughter(0).vz - daughter(1).vz"),
        probeMultiplicity = cms.InputTag("probeMultiplicity"),
    ),
    pairFlags = cms.PSet(),
    isMC           = cms.bool(False),
    addRunLumiInfo = cms.bool(True),
)


process.load("MuonAnalysis.TagAndProbe.muon.tag_probe_muon_extraIso_cfi")

process.tnpSimpleSequence = cms.Sequence(
    process.tagMuons +
    process.oneTag     +
    process.probeMuons +
    process.tpPairs    +
    process.onePair    +
    process.muonDxyPVdzmin +
    process.nverticesModule +
    process.tagProbeSeparation +
    process.computeCorrectedIso + 
    process.probeMultiplicity + 
    process.splitTrackTagger +
    process.l1rate +
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

process.load("MuonAnalysis.TagAndProbe.tracking_reco_info_cff")

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
        Mu5_L2Mu3_Jpsi_L2 = LowPtTriggerFlagsEfficienciesProbe.Mu5_L2Mu3_Jpsi_L2,
        TM  = cms.string("isTrackerMuon"),
        Glb = cms.string("isGlobalMuon"),
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
    ),
    tagFlags = cms.PSet(
        Mu5_L2Mu3_Jpsi_MU = LowPtTriggerFlagsEfficienciesTag.Mu5_L2Mu3_Jpsi_MU,
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
    process.l1rate +
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

if True: # turn on for tracking efficiency from RECO/AOD + earlyGeneralTracks
    process.tracksNoMuonSeeded = cms.EDFilter("TrackSelector",
      src = cms.InputTag("generalTracks"),
      cut = cms.string(" || ".join("isAlgoInMask('%s')" % a for a in [
                    'initialStep', 'lowPtTripletStep', 'pixelPairStep', 'detachedTripletStep',
                    'mixedTripletStep', 'pixelLessStep', 'tobTecStep', 'jetCoreRegionalStep' ] ) )
    )
    process.pCutTracks0 = process.pCutTracks.clone(src = 'tracksNoMuonSeeded')
    process.tkTracks0 = process.tkTracks.clone(src = 'pCutTracks0')
    process.tkTracksNoJPsi0 = process.tkTracksNoJPsi.clone(src = 'tkTracks0')
    process.tkTracksNoBestJPsi0 = process.tkTracksNoBestJPsi.clone(src = 'tkTracks0')
    process.preTkMatchSequenceJPsi.replace(
            process.tkTracksNoJPsi, process.tkTracksNoJPsi +
            process.tracksNoMuonSeeded + process.pCutTracks0 + process.tkTracks0 + process.tkTracksNoJPsi0 +process.tkTracksNoBestJPsi0
    )
    process.staToTkMatch0 = process.staToTkMatch.clone(matched = 'tkTracks0')
    process.staToTkMatchNoJPsi0 = process.staToTkMatchNoJPsi.clone(matched = 'tkTracksNoJPsi0')
    process.staToTkMatchNoBestJPsi0 = process.staToTkMatchNoBestJPsi.clone(matched = 'tkTracksNoJPsi0')
    process.staToTkMatchSequenceJPsi.replace( process.staToTkMatch, process.staToTkMatch + process.staToTkMatch0 )
    process.staToTkMatchSequenceJPsi.replace( process.staToTkMatchNoJPsi, process.staToTkMatchNoJPsi + process.staToTkMatchNoJPsi0 )
    process.staToTkMatchSequenceJPsi.replace( process.staToTkMatchNoBestJPsi, process.staToTkMatchNoBestJPsi + process.staToTkMatchNoBestJPsi0 )
    process.tpTreeSta.variables.tk0_deltaR     = cms.InputTag("staToTkMatch0","deltaR")
    process.tpTreeSta.variables.tk0_deltaEta   = cms.InputTag("staToTkMatch0","deltaEta")
    process.tpTreeSta.variables.tk0_deltaR_NoJPsi   = cms.InputTag("staToTkMatchNoJPsi0","deltaR")
    process.tpTreeSta.variables.tk0_deltaEta_NoJPsi = cms.InputTag("staToTkMatchNoJPsi0","deltaEta")
    process.tpTreeSta.variables.tk0_deltaR_NoBestJPsi   = cms.InputTag("staToTkMatchNoBestJPsi0","deltaR")
    process.tpTreeSta.variables.tk0_deltaEta_NoBestJPsi = cms.InputTag("staToTkMatchNoBestJPsi0","deltaEta")

process.tagAndProbeSta = cms.Path( 
    process.fastFilter +
    process.HLTMu      +
    process.muonsSta                       +
    process.patMuonsWithTriggerSequenceSta +
    process.tnpSimpleSequenceSta
)



if False: # turn on for tracking efficiency using L1 seeds
    process.probeL1 = cms.EDFilter("CandViewSelector",
        src = cms.InputTag("l1extraParticles"),
        cut = cms.string("pt >= 2 && abs(eta) < 2.4"),
    )
    process.tpPairsTkL1 = process.tpPairs.clone(decay = "tagMuons@+ probeL1@-", cut = 'mass > 2')
    process.l1ToTkMatch    = process.staToTkMatch.clone(src = "probeL1", srcTrack="none")
    process.l1ToTkMatchNoJPsi = process.staToTkMatchNoJPsi.clone(src = "probeL1", srcTrack="none")
    process.l1ToTkMatchNoBestJPsi = process.staToTkMatchNoBestJPsi.clone(src = "probeL1", srcTrack="none")
    process.l1ToTkMatch0    = process.staToTkMatch0.clone(src = "probeL1", srcTrack="none")
    process.l1ToTkMatchNoJPsi0 = process.staToTkMatchNoJPsi0.clone(src = "probeL1", srcTrack="none")
    process.l1ToTkMatchNoBestJPsi0 = process.staToTkMatchNoBestJPsi0.clone(src = "probeL1", srcTrack="none")
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
            tk_deltaR_NoJPsi   = cms.InputTag("l1ToTkMatchNoJPsi","deltaR"),
            tk_deltaEta_NoJPsi = cms.InputTag("l1ToTkMatchNoJPsi","deltaEta"),
            tk_deltaR_NoBestJPsi   = cms.InputTag("l1ToTkMatchNoBestJPsi","deltaR"),
            tk_deltaEta_NoBestJPsi = cms.InputTag("l1ToTkMatchNoBestJPsi","deltaEta"),
            ## track matching variables (early general tracks)
            tk0_deltaR     = cms.InputTag("l1ToTkMatch0","deltaR"),
            tk0_deltaEta   = cms.InputTag("l1ToTkMatch0","deltaEta"),
            tk0_deltaR_NoJPsi   = cms.InputTag("l1ToTkMatchNoJPsi0","deltaR"),
            tk0_deltaEta_NoJPsi = cms.InputTag("l1ToTkMatchNoJPsi0","deltaEta"),
            tk0_deltaR_NoBestJPsi   = cms.InputTag("l1ToTkMatchNoBestJPsi0","deltaR"),
            tk0_deltaEta_NoBestJPsi = cms.InputTag("l1ToTkMatchNoBestJPsi0","deltaEta"),
        ),
        flags = cms.PSet(
        ),
        tagVariables = cms.PSet(
            pt = cms.string("pt"),
            eta = cms.string("eta"),
            phi = cms.string("phi"),
            nVertices   = cms.InputTag("nverticesModule"),
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
        process.preTkMatchSequenceJPsi +
        process.l1ToTkMatch + 
        process.l1ToTkMatchNoJPsi + process.l1ToTkMatchNoBestJPsi +
        process.l1ToTkMatch0 + 
        process.l1ToTkMatchNoJPsi0 + process.l1ToTkMatchNoBestJPsi0 +
        process.tpTreeL1
    )

process.schedule = cms.Schedule(
   process.tagAndProbe,
   process.tagAndProbeSta,
   #process.tagAndProbeTkL1
)

process.RandomNumberGeneratorService.tkTracksNoJPsi      = cms.PSet( initialSeed = cms.untracked.uint32(81) )
process.RandomNumberGeneratorService.tkTracksNoJPsi0      = cms.PSet( initialSeed = cms.untracked.uint32(81) )
process.RandomNumberGeneratorService.tkTracksNoBestJPsi  = cms.PSet( initialSeed = cms.untracked.uint32(81) )
process.RandomNumberGeneratorService.tkTracksNoBestJPsi0  = cms.PSet( initialSeed = cms.untracked.uint32(81) )

process.TFileService = cms.Service("TFileService", fileName = cms.string("tnpJPsi_Data.root"))

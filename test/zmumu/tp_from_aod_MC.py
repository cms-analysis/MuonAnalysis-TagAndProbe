import FWCore.ParameterSet.Config as cms

process = cms.Process("TagProbe")

process.load('FWCore.MessageService.MessageLogger_cfi')
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.source = cms.Source("PoolSource", 
    fileNames = cms.untracked.vstring(
	'/store/relval/CMSSW_3_8_0/RelValZMM/GEN-SIM-RECO/START38_V7-v1/0005/5AE92533-3D95-DF11-8784-003048678F9C.root',
	'/store/relval/CMSSW_3_8_0/RelValZMM/GEN-SIM-RECO/START38_V7-v1/0004/78DF3A0F-F794-DF11-8920-003048678B88.root',
	'/store/relval/CMSSW_3_8_0/RelValZMM/GEN-SIM-RECO/START38_V7-v1/0004/68666746-F994-DF11-B168-003048D15E24.root',
	'/store/relval/CMSSW_3_8_0/RelValZMM/GEN-SIM-RECO/START38_V7-v1/0004/4CC4F3F5-F894-DF11-8A07-003048678FAE.root',
	'/store/relval/CMSSW_3_8_0/RelValZMM/GEN-SIM-RECO/START38_V7-v1/0004/30A3A23A-F794-DF11-84F0-003048678D52.root',
	'/store/relval/CMSSW_3_8_0/RelValZMM/GEN-SIM-RECO/START38_V7-v1/0004/182179E1-1895-DF11-968A-0018F3D0970E.root',
    ),
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )    

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.GlobalTag.globaltag = cms.string('GR_R_38X_V8::All')

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
from MuonAnalysis.MuonAssociators.patMuonsWithTrigger_cff import *
changeRecoMuonInput(process, "mergedMuons")

VBTF_CUTS = ("numberOfMatches > 1 && muonID('GlobalMuonPromptTight') && abs(dB) < 0.2 && "+
             "track.numberOfValidHits > 10 && track.hitPattern.numberOfValidPixelHits > 0")
process.tagMuons = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("patMuonsWithTrigger"),
    cut = cms.string("pt > 15 && "+VBTF_CUTS+" && !triggerObjectMatchesByFilter('hltSingleMu9L3Filtered9').empty()"),
)

process.probeMuons = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("patMuonsWithTrigger"),
    cut = cms.string("track.isNonnull"),
)

process.tpPairs = cms.EDProducer("CandViewShallowCloneCombiner",
    cut = cms.string('40 < mass < 140'),
    decay = cms.string('tagMuons@+ probeMuons@-')
)

process.tagMuonsMCMatch = cms.EDProducer("MCTruthDeltaRMatcherNew",
    src = cms.InputTag("tagMuons"),
    matched = cms.InputTag("genParticles"),
    pdgId = cms.vint32(13),
    distMin = cms.double(0.3),
)
process.probeMuonsMCMatch = process.tagMuonsMCMatch.clone(src = "probeMuons")

process.tpTree = cms.EDAnalyzer("TagProbeFitTreeProducer",
    # choice of tag and probe pairs, and arbitration
    tagProbePairs = cms.InputTag("tpPairs"),
    arbitration   = cms.string("OneProbe"),
    # probe variables 
    variables = cms.PSet(
        pt  = cms.string("pt"),
        p   = cms.string("p"),
        eta = cms.string("eta"),
        phi = cms.string("phi"),
        abseta = cms.string("abs(eta)"),
        tkIso  = cms.string("isolationR03.sumPt"),
        combRelIso = cms.string("(isolationR03.emEt + isolationR03.hadEt + isolationR03.sumPt)/pt"),
    ),
    flags = cms.PSet(
       ## MuonID cuts
       Calo   = cms.string("isCaloMuon"),
       Glb    = cms.string("isGlobalMuon"),
       GlbPT  = cms.string("muonID('GlobalMuonPromptTight')"),
       TMA    = cms.string("muonID('TrackerMuonArbitrated')"),
       TMLSAT = cms.string("muonID('TMLastStationAngTight')"),
       VBTF   = cms.string(VBTF_CUTS),
       ## Isolation
       Isol    = cms.string("(isolationR03.emEt + isolationR03.hadEt + isolationR03.sumPt)/pt < 0.15"), 
       IsolTk3 = cms.string("isolationR03.sumPt < 3"), 
       ## Trigger
       Mu9       = cms.string("!triggerObjectMatchesByFilter('hltSingleMu9L3Filtered9').empty()"),
       DoubleMu3 = cms.string("(!triggerObjectMatchesByFilter('hltDiMuonL3PreFiltered').empty())"),
    ),
    isMC           = cms.bool(True),
    tagMatches       = cms.InputTag("tagMuonsMCMatch"),
    probeMatches     = cms.InputTag("probeMuonsMCMatch"),
    motherPdgId      = cms.vint32(22, 23),
    makeMCUnbiasTree       = cms.bool(True),
    checkMotherInUnbiasEff = cms.bool(True),
    allProbes              = cms.InputTag("probeMuons"),
)
process.tnpSimpleSequence = cms.Sequence(
    process.tagMuons   * process.tagMuonsMCMatch   *
    process.probeMuons * process.probeMuonsMCMatch *
    process.tpPairs    *
    process.tpTree
)

process.tagAndProbe = cms.Path( 
    process.mergedMuons                 *
    process.patMuonsWithTriggerSequence *
    process.tnpSimpleSequence
)

process.TFileService = cms.Service("TFileService", fileName = cms.string("tnpZ_MC.root"))

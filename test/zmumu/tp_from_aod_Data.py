import FWCore.ParameterSet.Config as cms

process = cms.Process("TagProbe")

process.load('FWCore.MessageService.MessageLogger_cfi')
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.source = cms.Source("PoolSource", 
    fileNames = cms.untracked.vstring(
	'file:/nfs/bluearc/group/skims/ewk/data/mu/Muv4/ewkMuSkim.Muv4_1_1_KZ9.root',
	'file:/nfs/bluearc/group/skims/ewk/data/mu/Muv4/ewkMuSkim.Muv4_2_1_7Ou.root',
	'file:/nfs/bluearc/group/skims/ewk/data/mu/Muv4/ewkMuSkim.Muv4_3_1_J0y.root',
	'file:/nfs/bluearc/group/skims/ewk/data/mu/Muv4/ewkMuSkim.Muv4_4_1_Kzq.root',
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
    isMC           = cms.bool(False),
    addRunLumiInfo = cms.bool(True),
)
process.tnpSimpleSequence = cms.Sequence(
    process.tagMuons   *
    process.probeMuons *
    process.tpPairs    *
    process.tpTree
)

from HLTrigger.HLTfilters.hltHighLevelDev_cfi import hltHighLevelDev
process.fastFilter = hltHighLevelDev.clone(HLTPaths = ['HLT_Mu9'], HLTPathsPrescales = [1])
process.tagAndProbe = cms.Path( 
    process.fastFilter                  +
    process.mergedMuons                 *
    process.patMuonsWithTriggerSequence *
    process.tnpSimpleSequence
)

process.TFileService = cms.Service("TFileService", fileName = cms.string("tnpZ.root"))

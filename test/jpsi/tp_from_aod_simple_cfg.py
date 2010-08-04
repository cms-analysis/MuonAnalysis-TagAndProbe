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

## Here we use Mu3_Track0_Jpsi
TAG_TRIGGER = ("(!triggerObjectMatchesByCollection('hltL3MuonCandidates').empty() && "+
               "  triggerObjectMatchesByCollection('hltL3MuonCandidates').at(0).hasFilterLabel('hltMu3TrackJpsiTrackMassFiltered'))")
process.tagMuons = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("patMuonsWithTrigger"),
    cut = cms.string("isGlobalMuon && "+TAG_TRIGGER),
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
    # probe variables 
    variables = cms.PSet(
        pt  = cms.string("pt"),
        p   = cms.string("p"),
        eta = cms.string("eta"),
        phi = cms.string("phi"),
        abseta = cms.string("abs(eta)"),
    ),
    flags = cms.PSet(
       ## Track quality cuts
       Tk_HP   = cms.string("track.quality('highPurity')"),
       Tk_QTF  = cms.string("track.numberOfValidHits > 11 && track.hitPattern.pixelLayersWithMeasurement > 1 && track.normalizedChi2 < 4 && abs(dB) < 3 && abs(track.dz) < 30"),
       Tk_VBTF = cms.string("track.numberOfValidHits > 10 && track.hitPattern.numberOfValidPixelHits > 0"),
       ## MuonID cuts
       Calo   = cms.string("isCaloMuon"),
       Glb    = cms.string("isGlobalMuon"),
       GlbPT  = cms.string("muonID('GlobalMuonPromptTight')"),
       TMA    = cms.string("muonID('TrackerMuonArbitrated')"),
       TMLSAT = cms.string("muonID('TMLastStationAngTight')"),
       VBTF   = cms.string("numberOfMatches > 1 && muonID('GlobalMuonPromptTight') && abs(dB) < 0.2"),
       ## Trigger (just a few examples)
       L1DoubleMuOpen       = cms.string("!triggerObjectMatchesByFilter('hltDoubleMuLevel1PathL1OpenFiltered').empty()"),
       L1DoubleMuOpen_Tight = cms.string("!triggerObjectMatchesByFilter('hltL1DoubleMuOpenTightL1Filtered').empty()"),
       DoubleMu0          = cms.string("!triggerObjectMatchesByFilter('hltDiMuonL3PreFiltered0').empty()"),
       DoubleMu3          = cms.string("!triggerObjectMatchesByFilter('hltDiMuonL3PreFiltered').empty()"),
       Mu0_Track0_Jpsi_MU = cms.string("!triggerObjectMatchesByCollection('hltL3MuonCandidates').empty() && "+
                                       " triggerObjectMatchesByCollection('hltL3MuonCandidates').at(0).hasFilterLabel('hltMu0TrackJpsiTrackMassFiltered')"),
       Mu0_Track0_Jpsi_TK = cms.string("!triggerObjectMatchesByCollection('hltMuTrackJpsiCtfTrackCands').empty() && "+
                                       " triggerObjectMatchesByCollection('hltMuTrackJpsiCtfTrackCands').at(0).hasFilterLabel('hltMu0TrackJpsiTrackMassFiltered')"),
       Mu3_Track0_Jpsi_MU = cms.string("!triggerObjectMatchesByCollection('hltL3MuonCandidates').empty() && "+
                                       " triggerObjectMatchesByCollection('hltL3MuonCandidates').at(0).hasFilterLabel('hltMu3TrackJpsiTrackMassFiltered')"),
       Mu3_Track0_Jpsi_TK = cms.string("!triggerObjectMatchesByCollection('hltMuTrackJpsiCtfTrackCands').empty() && "+
                                       " triggerObjectMatchesByCollection('hltMuTrackJpsiCtfTrackCands').at(0).hasFilterLabel('hltMu3TrackJpsiTrackMassFiltered')"),
       Mu5_Track0_Jpsi_MU = cms.string("!triggerObjectMatchesByCollection('hltL3MuonCandidates').empty() && "+
                                       " triggerObjectMatchesByCollection('hltL3MuonCandidates').at(0).hasFilterLabel('hltMu5TrackJpsiTrackMassFiltered')"),
       Mu5_Track0_Jpsi_TK = cms.string("!triggerObjectMatchesByCollection('hltMuTrackJpsiCtfTrackCands').empty() && "+
                                       " triggerObjectMatchesByCollection('hltMuTrackJpsiCtfTrackCands').at(0).hasFilterLabel('hltMu5TrackJpsiTrackMassFiltered')"),
       Mu0_TkMu0_Jpsi_MU = cms.string("!triggerObjectMatchesByCollection('hltL3MuonCandidates').empty() && "+
                                      " triggerObjectMatchesByCollection('hltL3MuonCandidates').at(0).hasFilterLabel('hltMu0TkMuJpsiTkMuMassFiltered')"),
       Mu0_TkMu0_Jpsi_TM = cms.string("!triggerObjectMatchesByCollection('hltMuTkMuJpsiTrackerMuonCands').empty() && "+
                                      " triggerObjectMatchesByCollection('hltMuTkMuJpsiTrackerMuonCands').at(0).hasFilterLabel('hltMu0TkMuJpsiTkMuMassFiltered')"),
       ## Acceptance definition
       Acc_JPsi = cms.string("(abs(eta) <= 1.3 && pt > 3.3) || (1.3 < abs(eta) <= 2.2 && p > 2.9) || (2.2 < abs(eta) <= 2.4  && pt > 0.8)"),
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
process.fastFilter = hltHighLevelDev.clone(HLTPaths = ['HLT_Mu3_Track0_Jpsi'], HLTPathsPrescales = [1])
process.tagAndProbe = cms.Path( 
    process.fastFilter                  +
    process.mergedMuons                 *
    process.patMuonsWithTriggerSequence *
    process.tnpSimpleSequence
)

process.TFileService = cms.Service("TFileService", fileName = cms.string("tnpSimple.root"))

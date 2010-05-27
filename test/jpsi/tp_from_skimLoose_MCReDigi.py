import FWCore.ParameterSet.Config as cms

process = cms.Process("TagProbe")

process.load('FWCore.MessageService.MessageLogger_cfi')
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.source = cms.Source("PoolSource", 
    fileNames = cms.untracked.vstring(
        'root://pcmssd12.cern.ch//data/gpetrucc/7TeV/tnp/ppMuX_Spring10ReDigi_skimJPsiLoose_v2/skimJPsiLoose_Spring10ReDigi_1_1.root',
        'root://pcmssd12.cern.ch//data/gpetrucc/7TeV/tnp/ppMuX_Spring10ReDigi_skimJPsiLoose_v2/skimJPsiLoose_Spring10ReDigi_2_1.root',
        'root://pcmssd12.cern.ch//data/gpetrucc/7TeV/tnp/ppMuX_Spring10ReDigi_skimJPsiLoose_v2/skimJPsiLoose_Spring10ReDigi_3_1.root',
        'root://pcmssd12.cern.ch//data/gpetrucc/7TeV/tnp/ppMuX_Spring10ReDigi_skimJPsiLoose_v2/skimJPsiLoose_Spring10ReDigi_4_1.root',
        'root://pcmssd12.cern.ch//data/gpetrucc/7TeV/tnp/ppMuX_Spring10ReDigi_skimJPsiLoose_v2/skimJPsiLoose_Spring10ReDigi_5_2.root',
    )
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )    

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Geometry.CommonDetUnit.globalTrackingGeometry_cfi")
process.load("TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorAny_cfi")
process.load("TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorAlong_cfi")
process.load("TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorOpposite_cfi")
process.GlobalTag.globaltag = cms.string('MC_3XY_V25::All')

process.load("MuonAnalysis.TagAndProbe.jpsi.tp_from_skim_Tracking_cff")
process.load("MuonAnalysis.TagAndProbe.jpsi.tp_from_skim_MuonID_cff")
process.load("MuonAnalysis.TagAndProbe.jpsi.tp_from_skim_Trigger_cff")

from MuonAnalysis.TagAndProbe.jpsi.tp_from_skim_common_cff import *

## Loosen tag requirements to the minimum
process.tagMuons1Mu.cut = "isGlobalMuon && !triggerObjectMatchesByFilter('hltSingleMu3L3Filtered3').empty() && " + TRACK_CUTS;
process.tagMuons2Mu.cut = "isGlobalMuon && !triggerObjectMatchesByFilter('hltSingleMu3L3Filtered3').empty() && " + TRACK_CUTS;
#process.tagMuons2Mu.cut = "isGlobalMuon && !triggerObjectMatchesByFilter('hltDoubleMuLevel1PathL1OpenFiltered').empty() && " + TRACK_CUTS;

from MuonAnalysis.TagAndProbe.jpsi.tp_from_skim_common_cff import addDiMuonSeparationVariables
addDiMuonSeparationVariables(process, process.tnpSequenceTrigger, process.histoTrigger)
addDiMuonSeparationVariables(process, process.tnpSequenceMuonID,  process.histoMuFromTk)
addDiMuonSeparationVariables(process, process.tnpSequenceMuonID,  process.histoMuFromCal)

## On REDIGI
for K,V in allTPTreeProducers(process): 
    V.tagFlags.Mu0_Track0_JPsi = cms.string("!triggerObjectMatchesByFilter('hltMu0TrackJpsiTrackMassFiltered').empty() && !triggerObjectMatchesByCollection('hltL3MuonCandidates::REDIGI').empty()")
    V.tagFlags.Mu3_Track0_JPsi = cms.string("!triggerObjectMatchesByFilter('hltMu3TrackJpsiTrackMassFiltered').empty() && !triggerObjectMatchesByCollection('hltL3MuonCandidates::REDIGI').empty()")

process.tagAndProbe = cms.Path( 
    process.tnpCommonSequence    *
    ( process.tnpSequenceTracking  +
      process.tnpSequenceMuonID    +
      process.tnpSequenceTrigger   )
)


##     ___        _               _     _   _ _     _            
##    / _ \ _   _| |_ _ __  _   _| |_  | | | (_)___| |_ ___  ___ 
##   | | | | | | | __| '_ \| | | | __| | |_| | / __| __/ _ \/ __|
##   | |_| | |_| | |_| |_) | |_| | |_  |  _  | \__ \ || (_) \__ \
##    \___/ \__,_|\__| .__/ \__,_|\__| |_| |_|_|___/\__\___/|___/
##                   |_|                                         
##   
process.TFileService = cms.Service("TFileService", fileName = cms.string("tnpJPsi_MC.root"))



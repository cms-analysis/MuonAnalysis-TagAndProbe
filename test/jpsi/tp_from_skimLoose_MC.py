import FWCore.ParameterSet.Config as cms

process = cms.Process("TagProbe")

process.load('FWCore.MessageService.MessageLogger_cfi')
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.source = cms.Source("PoolSource", 
    fileNames = cms.untracked.vstring(
        'rfio:/castor/cern.ch/user/g/gpetrucc/7TeV/MC/TnP/JPsiToMuMu_Spring10_skimJPsiLoose_v1/skimJPsiLoose_ppxMuXLoose900GeV_1_1.root',
        'rfio:/castor/cern.ch/user/g/gpetrucc/7TeV/MC/TnP/JPsiToMuMu_Spring10_skimJPsiLoose_v1/skimJPsiLoose_ppxMuXLoose900GeV_2_0.root',
        'rfio:/castor/cern.ch/user/g/gpetrucc/7TeV/MC/TnP/JPsiToMuMu_Spring10_skimJPsiLoose_v1/skimJPsiLoose_ppxMuXLoose900GeV_3_0.root',
        'rfio:/castor/cern.ch/user/g/gpetrucc/7TeV/MC/TnP/JPsiToMuMu_Spring10_skimJPsiLoose_v1/skimJPsiLoose_ppxMuXLoose900GeV_4_0.root',
        #'root://pcmssd12.cern.ch//data/gpetrucc/7TeV/tnp/JPsiToMuMu_Spring10_skimJPsiLoose_v1/skimJPsiLoose_ppxMuXLoose900GeV_1_1.root',
        #'root://pcmssd12.cern.ch//data/gpetrucc/7TeV/tnp/JPsiToMuMu_Spring10_skimJPsiLoose_v1/skimJPsiLoose_ppxMuXLoose900GeV_2_0.root',
        #'root://pcmssd12.cern.ch//data/gpetrucc/7TeV/tnp/JPsiToMuMu_Spring10_skimJPsiLoose_v1/skimJPsiLoose_ppxMuXLoose900GeV_3_0.root',
        #'root://pcmssd12.cern.ch//data/gpetrucc/7TeV/tnp/JPsiToMuMu_Spring10_skimJPsiLoose_v1/skimJPsiLoose_ppxMuXLoose900GeV_4_0.root',
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
process.tagMuons1Mu.cut = PASSING_GLB_CUT + " && !triggerObjectMatchesByFilter('hltL1MuOpenL1Filtered0').empty() && " + TRACK_CUTS;
process.tagMuons2Mu.cut = PASSING_GLB_CUT + " && !triggerObjectMatchesByFilter('hltL1MuOpenL1Filtered0').empty() && " + TRACK_CUTS;

## On REDIGI
#for K,V in allTPTreeProducers(process): 
#    V.tagFlags.HLTMu0Tk = cms.string("!triggerObjectMatchesByFilter('hltMu0TrackJpsiTrackMassFiltered').empty() && !triggerObjectMatchesByCollection('hltL3MuonCandidates::REDIGI').empty()")
#    V.tagFlags.HLTMu3Tk = cms.string("!triggerObjectMatchesByFilter('hltMu3TrackJpsiTrackMassFiltered').empty() && !triggerObjectMatchesByCollection('hltL3MuonCandidates::REDIGI').empty()")

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



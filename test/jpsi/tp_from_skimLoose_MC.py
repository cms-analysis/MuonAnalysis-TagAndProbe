import FWCore.ParameterSet.Config as cms

process = cms.Process("TagProbe")

process.load('FWCore.MessageService.MessageLogger_cfi')
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.source = cms.Source("PoolSource", 
    fileNames = cms.untracked.vstring(
        'root://pcmssd12.cern.ch//data/gpetrucc/7TeV/tnp/JPsiToMuMu_Spring10_skimJPsiLoose_v3/skimJPsiLoose_Spring10_1_1.root',
        'root://pcmssd12.cern.ch//data/gpetrucc/7TeV/tnp/JPsiToMuMu_Spring10_skimJPsiLoose_v3/skimJPsiLoose_Spring10_2_1.root',
        'root://pcmssd12.cern.ch//data/gpetrucc/7TeV/tnp/JPsiToMuMu_Spring10_skimJPsiLoose_v3/skimJPsiLoose_Spring10_3_1.root',
        'root://pcmssd12.cern.ch//data/gpetrucc/7TeV/tnp/JPsiToMuMu_Spring10_skimJPsiLoose_v3/skimJPsiLoose_Spring10_4_1.root',
        'root://pcmssd12.cern.ch//data/gpetrucc/7TeV/tnp/JPsiToMuMu_Spring10_skimJPsiLoose_v3/skimJPsiLoose_Spring10_5_1.root',
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

from MuonAnalysis.TagAndProbe.jpsi.tp_from_skim_common_cff import addDiMuonSeparationVariables
addDiMuonSeparationVariables(process, process.tnpSequenceTrigger, process.histoTrigger)
addDiMuonSeparationVariables(process, process.tnpSequenceMuonID,  process.histoMuFromTk)
addDiMuonSeparationVariables(process, process.tnpSequenceMuonID,  process.histoMuFromCal)

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

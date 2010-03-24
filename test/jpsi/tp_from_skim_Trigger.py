import FWCore.ParameterSet.Config as cms

process = cms.Process("TagProbe")

process.load('FWCore.MessageService.MessageLogger_cfi')
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.source = cms.Source("PoolSource", 
    fileNames = cms.untracked.vstring(
        'root://pcmssd12.cern.ch//data/gpetrucc/tnp/ppMuX_skimJPsi_1.root',
        'root://pcmssd12.cern.ch//data/gpetrucc/tnp/ppMuX_skimJPsi_2.root',
        'root://pcmssd12.cern.ch//data/gpetrucc/tnp/JPsiMuMu_skimJPsi_1.root',
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

process.load("MuonAnalysis.TagAndProbe.jpsi.tp_from_skim_Trigger_cff")

process.tagAndProbe = cms.Path(
    process.tnpCommonSequence *
    process.tnpSequenceTrigger
)


process.TFileService = cms.Service("TFileService", fileName = cms.string("tnpJPsi_Trigger.root"))

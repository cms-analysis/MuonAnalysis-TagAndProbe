import FWCore.ParameterSet.Config as cms

process = cms.Process("myprocess")
process.load("CondCore.DBCommon.CondDBCommon_cfi")

process.CondDBCommon.connect = 'sqlite_file:MuonPhysicsPerformance.db'

process.load ("MuonAnalysis.TagAndProbe.MuonPerformanceESProducer_cfi")
#
# change inside the source
#

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
    '/store/relval/CMSSW_3_1_2/RelValZMM/GEN-SIM-RECO/STARTUP31X_V2-v1/0007/EECE7AB6-CC78-DE11-805C-0019B9F709A4.root'    
#    '/store/relval/CMSSW_3_3_0_pre2/RelValZMM/GEN-SIM-RECO/STARTUP31X_V7-v1/0003/16950730-469C-DE11-9130-001731A281B1.root'        
    )
                            )


process.PoolDBESSource = cms.ESSource("PoolDBESSource",
                                      process.CondDBCommon,
                                      toGet = cms.VPSet(
    cms.PSet(
    record = cms.string('PerformanceWPRecord'),
    tag = cms.string('GLBMU_WP'),
    label = cms.untracked.string('GLBMU_WP')
    ),
    cms.PSet(
    record = cms.string('PerformancePayloadRecord'),
    tag = cms.string('GLBMU_TABLE'),
    label = cms.untracked.string('GLBMU_TABLE')
    )))

process.PoolDBESSource2 = cms.ESSource("PoolDBESSource",
                                      process.CondDBCommon,
                                      toGet = cms.VPSet(
    cms.PSet(
    record = cms.string('PerformanceWPRecord'),
    tag = cms.string('TRKEFFMU_WP'),
    label = cms.untracked.string('TRKEFFMU_WP')
    ),
    cms.PSet(
    record = cms.string('PerformancePayloadRecord'),
    tag = cms.string('TRKEFFMU_TABLE'),
    label = cms.untracked.string('TRKEFFMU_TABLE')
    )))

process.demo2 = cms.EDAnalyzer('MuTestPerformanceFW_ES',
                               AlgoName1 = cms.string('GlobalMuonFromTrackerTrackJpsi'),
                               AlgoName2 = cms.string('TrackerTrackFromStandaloneMuonJpsi'))

#
# change inside the source
#
process.MuonPerformanceESProducer_GlobalMuon.PayloadName = "GLBMU_TABLE"
process.MuonPerformanceESProducer_GlobalMuon.WorkingPointName = "GLBMU_WP"
process.MuonPerformanceESProducer_GlobalMuon.ComponentName = "GlobalMuonFromTrackerTrackJpsi"
process.MuonPerformanceESProducer_TrackerTrackMuon.PayloadName = "TRKEFFMU_TABLE"
process.MuonPerformanceESProducer_TrackerTrackMuon.WorkingPointName = "TRKEFFMU_WP"
process.MuonPerformanceESProducer_TrackerTrackMuon.ComponentName = "TrackerTrackFromStandaloneMuonJpsi"

process.p = cms.Path(process.demo2)

#print process.dumpPython()

#


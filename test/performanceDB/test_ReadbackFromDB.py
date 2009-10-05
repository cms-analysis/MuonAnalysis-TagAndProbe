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
    tag = cms.string('GLBMUZ_WP'),
    label = cms.untracked.string('GLBMUZ_WP')
    ),
    cms.PSet(
    record = cms.string('PerformancePayloadRecord'),
    tag = cms.string('GLBMUZ_TABLE'),
    label = cms.untracked.string('GLBMUZ_TABLE')
    )))

process.PoolDBESSource2 = cms.ESSource("PoolDBESSource",
                                      process.CondDBCommon,
                                      toGet = cms.VPSet(
    cms.PSet(
    record = cms.string('PerformanceWPRecord'),
    tag = cms.string('TRKEFFMUZ_WP'),
    label = cms.untracked.string('TRKEFFMUZ_WP')
    ),
    cms.PSet(
    record = cms.string('PerformancePayloadRecord'),
    tag = cms.string('TRKEFFMUZ_TABLE'),
    label = cms.untracked.string('TRKEFFMUZ_TABLE')
    )))

process.PoolDBESSource3 = cms.ESSource("PoolDBESSource",
                                       process.CondDBCommon,
                                       toGet = cms.VPSet(
    cms.PSet(
    record = cms.string('PerformanceWPRecord'),
    tag = cms.string('TRGMUZ_WP'),
    label = cms.untracked.string('TRGMUZ_WP')
    ),
    cms.PSet(
    record = cms.string('PerformancePayloadRecord'),
    tag = cms.string('TRGMUZ_TABLE'),
    label = cms.untracked.string('TRGMUZ_TABLE')
    )))

process.demo2 = cms.EDAnalyzer('MuTestPerformanceFW_ES',
                               AlgoNames = cms.vstring(
    'TriggerMuonFromGlobalMuonZ',
    'GlobalMuonFromTrackerTrackZ',
    'TrackerTrackFromStandaloneMuonZ'))

#
# change inside the source
#
process.MuonPerformanceESProducer_GlobalMuon.PayloadName = "GLBMUZ_TABLE"
process.MuonPerformanceESProducer_GlobalMuon.WorkingPointName = "GLBMUZ_WP"
process.MuonPerformanceESProducer_GlobalMuon.ComponentName = "GlobalMuonFromTrackerTrackZ"
process.MuonPerformanceESProducer_TrackerTrackMuon.PayloadName = "TRKEFFMUZ_TABLE"
process.MuonPerformanceESProducer_TrackerTrackMuon.WorkingPointName = "TRKEFFMUZ_WP"
process.MuonPerformanceESProducer_TrackerTrackMuon.ComponentName = "TrackerTrackFromStandaloneMuonZ"
process.MuonPerformanceESProducer_TriggerMuon.PayloadName = "TRGMUZ_TABLE"
process.MuonPerformanceESProducer_TriggerMuon.WorkingPointName = "TRGMUZ_WP"
process.MuonPerformanceESProducer_TriggerMuon.ComponentName = "TriggerMuonFromGlobalMuonZ"


process.p = cms.Path(process.demo2)

#print process.dumpPython()

#


import FWCore.ParameterSet.Config as cms

process = cms.Process("myprocess")
process.load("CondCore.DBCommon.CondDBCommon_cfi")

process.CondDBCommon.connect = 'sqlite_file:MuonPhysicsPerformance.db'

process.load ("MuonAnalysis.TagAndProbe.MuonPerformanceESProducer_cfi")
#
# change inside the source
#

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
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
    tag = cms.string('GLBMUJPSI_OCTXTEST_WP'),
    label = cms.untracked.string('GLBMUJPSI_OCTXTEST_WP')
    ),
    cms.PSet(
    record = cms.string('PerformancePayloadRecord'),
    tag = cms.string('GLBMUJPSI_OCTXTEST_TABLE'),
    label = cms.untracked.string('GLBMUJPSI_OCTXTEST_TABLE')
    ),
    cms.PSet(
    record = cms.string('PerformanceWPRecord'),
    tag = cms.string('GLBMUZ_OCTXTEST_WP'),
    label = cms.untracked.string('GLBMUZ_OCTXTEST_WP')
    ),
    cms.PSet(
    record = cms.string('PerformancePayloadRecord'),
    tag = cms.string('GLBMUZ_OCTXTEST_TABLE'),
    label = cms.untracked.string('GLBMUZ_OCTXTEST_TABLE')
    )))


process.PoolDBESSource2 = cms.ESSource("PoolDBESSource",
                                      process.CondDBCommon,
                                      toGet = cms.VPSet(
    cms.PSet(
    record = cms.string('PerformanceWPRecord'),
    tag = cms.string('TRKEFFMUJPSI_OCTXTEST_WP'),
    label = cms.untracked.string('TRKEFFMUJPSI_OCTXTEST_WP')
    ),
    cms.PSet(
    record = cms.string('PerformancePayloadRecord'),
    tag = cms.string('TRKEFFMUJPSI_OCTXTEST_TABLE'),
    label = cms.untracked.string('TRKEFFMUJPSI_OCTXTEST_TABLE')
    ),
    cms.PSet(
    record = cms.string('PerformanceWPRecord'),
    tag = cms.string('TRKEFFMUZ_OCTXTEST_WP'),
    label = cms.untracked.string('TRKEFFMUZ_OCTXTEST_WP')
    ),
    cms.PSet(
    record = cms.string('PerformancePayloadRecord'),
    tag = cms.string('TRKEFFMUZ_OCTXTEST_TABLE'),
    label = cms.untracked.string('TRKEFFMUZ_OCTXTEST_TABLE')
    )))

process.PoolDBESSource3 = cms.ESSource("PoolDBESSource",
                                       process.CondDBCommon,
                                       toGet = cms.VPSet(
    cms.PSet(
    record = cms.string('PerformanceWPRecord'),
    tag = cms.string('TRGMUJPSI_OCTXTEST_WP'),
    label = cms.untracked.string('TRGMUJPSI_OCTXTEST_WP')
    ),
    cms.PSet(
    record = cms.string('PerformancePayloadRecord'),
    tag = cms.string('TRGMUJPSI_OCTXTEST_TABLE'),
    label = cms.untracked.string('TRGMUJPSI_OCTXTEST_TABLE')
    ),
    cms.PSet(
    record = cms.string('PerformanceWPRecord'),
    tag = cms.string('TRGMUZ_OCTXTEST_WP'),
    label = cms.untracked.string('TRGMUZ_OCTXTEST_WP')
    ),
    cms.PSet(
    record = cms.string('PerformancePayloadRecord'),
    tag = cms.string('TRGMUZ_OCTXTEST_TABLE'),
    label = cms.untracked.string('TRGMUZ_OCTXTEST_TABLE')
    )))

# PAG-specific selections
process.PoolDBESSource4 =  cms.ESSource("PoolDBESSource",
                                       process.CondDBCommon,
                                       toGet = cms.VPSet(
    cms.PSet( 
    record = cms.string('PerformancePayloadRecord'), 
    tag = cms.string('TRGMUJPSI_JPSIANAL_OCTXTEST_TABLE'), 
    label = cms.untracked.string('TRGMUJPSI_JPSIANAL_OCTXTEST_TABLE') 
    ), 
    cms.PSet( 
    record = cms.string('PerformanceWPRecord'), 
    tag = cms.string('TRGMUJPSI_JPSIANAL_OCTXTEST_WP'), 
    label = cms.untracked.string('TRGMUJPSI_JPSIANAL_OCTXTEST_WP') 
    ), 
    cms.PSet( 
    record = cms.string('PerformancePayloadRecord'), 
    tag = cms.string('MUJPSI_JPSIGLBANAL_OCTXTEST_TABLE'), 
    label = cms.untracked.string('MUJPSI_JPSIGLBANAL_OCTXTEST_TABLE') 
    ), 
    cms.PSet( 
    record = cms.string('PerformanceWPRecord'), 
    tag = cms.string('MUJPSI_JPSIGLBANAL_OCTXTEST_WP'), 
    label = cms.untracked.string('MUJPSI_JPSIGLBANAL_OCTXTEST_WP') 
    ),     
    cms.PSet( 
    record = cms.string('PerformancePayloadRecord'), 
    tag = cms.string('MUJPSI_JPSITKMANAL_OCTXTEST_TABLE'), 
    label = cms.untracked.string('MUJPSI_JPSITKMANAL_OCTXTEST_TABLE') 
    ), 
    cms.PSet( 
    record = cms.string('PerformanceWPRecord'), 
    tag = cms.string('MUJPSI_JPSITKMANAL_OCTXTEST_WP'), 
    label = cms.untracked.string('MUJPSI_JPSITKMANAL_OCTXTEST_WP') 
    ), 
    cms.PSet( 
    record = cms.string('PerformancePayloadRecord'), 
    tag = cms.string('MUJPSI_JPSIPLUSMUANAL_OCTXTEST_TABLE'), 
    label = cms.untracked.string('MUJPSI_JPSIPLUSMUANAL_OCTXTEST_TABLE') 
    ), 
    cms.PSet( 
    record = cms.string('PerformanceWPRecord'), 
    tag = cms.string('MUJPSI_JPSIPLUSMUANAL_OCTXTEST_WP'), 
    label = cms.untracked.string('MUJPSI_JPSIPLUSMUANAL_OCTXTEST_WP') 
    ),     
    cms.PSet( 
    record = cms.string('PerformancePayloadRecord'), 
    tag = cms.string('MUJPSI_BEXCLANAL_OCTXTEST_TABLE'), 
    label = cms.untracked.string('MUJPSI_BEXCLANAL_OCTXTEST_TABLE') 
    ), 
    cms.PSet( 
    record = cms.string('PerformanceWPRecord'), 
    tag = cms.string('MUJPSI_BEXCLANAL_OCTXTEST_WP'), 
    label = cms.untracked.string('MUJPSI_BEXCLANAL_OCTXTEST_WP') 
    )))


#
# change inside the source
#
process.MuonPerformanceESProducer_GlobalMuon1.PayloadName = "GLBMUZ_OCTXTEST_TABLE"
process.MuonPerformanceESProducer_GlobalMuon1.WorkingPointName = "GLBMUZ_OCTXTEST_WP"
process.MuonPerformanceESProducer_GlobalMuon1.ComponentName = "GlobalMuonFromTrackerTrackZ"
process.MuonPerformanceESProducer_TrackerTrackMuon1.PayloadName = "TRKEFFMUZ_OCTXTEST_TABLE"
process.MuonPerformanceESProducer_TrackerTrackMuon1.WorkingPointName = "TRKEFFMUZ_OCTXTEST_WP"
process.MuonPerformanceESProducer_TrackerTrackMuon1.ComponentName = "TrackerTrackFromStandaloneMuonZ"
process.MuonPerformanceESProducer_TriggerMuon1.PayloadName = "TRGMUZ_OCTXTEST_TABLE"
process.MuonPerformanceESProducer_TriggerMuon1.WorkingPointName = "TRGMUZ_OCTXTEST_WP"
process.MuonPerformanceESProducer_TriggerMuon1.ComponentName = "TriggerMuonFromGlobalMuonZ"

process.MuonPerformanceESProducer_GlobalMuon2.PayloadName = "GLBMUJPSI_OCTXTEST_TABLE"
process.MuonPerformanceESProducer_GlobalMuon2.WorkingPointName = "GLBMUJPSI_OCTXTEST_WP"
process.MuonPerformanceESProducer_GlobalMuon2.ComponentName = "GlobalMuonFromTrackerTrackJpsi"
process.MuonPerformanceESProducer_TrackerTrackMuon2.PayloadName = "TRKEFFMUJPSI_OCTXTEST_TABLE"
process.MuonPerformanceESProducer_TrackerTrackMuon2.WorkingPointName = "TRKEFFMUJPSI_OCTXTEST_WP"
process.MuonPerformanceESProducer_TrackerTrackMuon2.ComponentName = "TrackerTrackFromStandaloneMuonJpsi"
process.MuonPerformanceESProducer_TriggerMuon2.PayloadName = "TRGMUJPSI_OCTXTEST_TABLE"
process.MuonPerformanceESProducer_TriggerMuon2.WorkingPointName = "TRGMUJPSI_OCTXTEST_WP"
process.MuonPerformanceESProducer_TriggerMuon2.ComponentName = "TriggerMuonFromGlobalMuonJpsi"

# PAG-specific selections
process.MuonPerformanceESProducer_TriggerMuon2.PayloadName = "TRGMUJPSI_JPSIANAL_OCTXTEST_TABLE"
process.MuonPerformanceESProducer_TriggerMuon2.WorkingPointName = "TRGMUJPSI_JPSIANAL_OCTXTEST_WP"
process.MuonPerformanceESProducer_TriggerMuon2.ComponentName = "TriggerMuonFromGlobalMuonJpsi_JpsiAnal"
process.MuonPerformanceESProducer_Muon1.PayloadName = "MUJPSI_JPSIGLBANAL_OCTXTEST_TABLE"
process.MuonPerformanceESProducer_Muon1.WorkingPointName = "MUJPSI_JPSIGLBANAL_OCTXTEST_WP"
process.MuonPerformanceESProducer_Muon1.ComponentName = "MuonFromTrackerTrackJpsi_JpsiGlbAnal"
process.MuonPerformanceESProducer_Muon2.PayloadName = "MUJPSI_JPSITKMANAL_OCTXTEST_TABLE"
process.MuonPerformanceESProducer_Muon2.WorkingPointName = "MUJPSI_JPSITKMANAL_OCTXTEST_WP"
process.MuonPerformanceESProducer_Muon2.ComponentName = "MuonFromTrackerTrackJpsi_JpsiTkMAnal"
process.MuonPerformanceESProducer_Muon3.PayloadName = "MUJPSI_JPSIPLUSMUANAL_OCTXTEST_TABLE"
process.MuonPerformanceESProducer_Muon3.WorkingPointName = "MUJPSI_JPSIPLUSMUANAL_OCTXTEST_WP"
process.MuonPerformanceESProducer_Muon3.ComponentName = "MuonFromTrackerTrackJpsi_JpsiPlusMuAnal"
process.MuonPerformanceESProducer_Muon3.PayloadName = "MUJPSI_BEXCLANAL_OCTXTEST_TABLE"
process.MuonPerformanceESProducer_Muon3.WorkingPointName = "MUJPSI_BEXCLANAL_OCTXTEST_WP"
process.MuonPerformanceESProducer_Muon3.ComponentName = "MuonFromTrackerTrackJpsi_BExclAnal"


process.demo2 = cms.EDAnalyzer('MuTestPerformanceFW_ES',
                               AlgoNames = cms.vstring(
    'TriggerMuonFromGlobalMuonZ',
    'GlobalMuonFromTrackerTrackZ',
    'TrackerTrackFromStandaloneMuonZ',
    'TriggerMuonFromGlobalMuonJpsi',
    'GlobalMuonFromTrackerTrackJpsi',
    'TrackerTrackFromStandaloneMuonJpsi',
    'TriggerMuonFromGlobalMuonJpsi_JpsiAnal', 
    'MuonFromTrackerTrackJpsi_JpsiGlbAnal', 
    'MuonFromTrackerTrackJpsi_JpsiTkMAnal', 
    'MuonFromTrackerTrackJpsi_JpsiPlusMuAnal', 
    'MuonFromTrackerTrackJpsi_BExclAnal'
    ))

process.p = cms.Path(process.demo2)

#print process.dumpPython()

#


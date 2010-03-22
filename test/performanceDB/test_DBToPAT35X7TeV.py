import FWCore.ParameterSet.Config as cms

process = cms.Process("MUONEFF")
process.load("CondCore.DBCommon.CondDBCommon_cfi")

process.CondDBCommon.connect = 'sqlite_file:MuonPhysicsPerformance7TeV.db'

process.load ("MuonAnalysis.TagAndProbe.MuonPerformanceESProducer_cfi")
#
# change inside the source
#

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
    '/store/relval/CMSSW_3_5_3/RelValJpsiMM/GEN-SIM-RECO/START3X_V24-v1/0002/D898733A-3C28-DF11-A76C-002618943979.root'
    )
                            )


process.PoolDBESSource = cms.ESSource("PoolDBESSource",
                                      process.CondDBCommon,
                                      toGet = cms.VPSet(
    cms.PSet(
    record = cms.string('PerformanceWPRecord'),
    tag = cms.string('GLBMUJPSI_TEST7TEV_WP'),
    label = cms.untracked.string('GLBMUJPSI_TEST7TEV_WP')
    ),
    cms.PSet(
    record = cms.string('PerformancePayloadRecord'),
    tag = cms.string('GLBMUJPSI_TEST7TEV_TABLE'),
    label = cms.untracked.string('GLBMUJPSI_TEST7TEV_TABLE')
    )))


# change inside the source

process.MuonPerformanceESProducer_GlobalMuon1.PayloadName = "GLBMUZ_TEST7TEV_TABLE"
process.MuonPerformanceESProducer_GlobalMuon1.WorkingPointName = "GLBMUZ_TEST7TEV_WP"
process.MuonPerformanceESProducer_GlobalMuon1.ComponentName = "GlobalMuonFromTrackerTrackZ"
process.MuonPerformanceESProducer_TrackerTrackMuon1.PayloadName = "TRKEFFMUZ_TEST7TEV_TABLE"
process.MuonPerformanceESProducer_TrackerTrackMuon1.WorkingPointName = "TRKEFFMUZ_TEST7TEV_WP"
process.MuonPerformanceESProducer_TrackerTrackMuon1.ComponentName = "TrackerTrackFromStandaloneMuonZ"
process.MuonPerformanceESProducer_TriggerMuon1.PayloadName = "TRGMUZ_TEST7TEV_TABLE"
process.MuonPerformanceESProducer_TriggerMuon1.WorkingPointName = "TRGMUZ_TEST7TEV_WP"
process.MuonPerformanceESProducer_TriggerMuon1.ComponentName = "TriggerMuonFromGlobalMuonZ"
process.MuonPerformanceESProducer_GlobalMuon2.PayloadName = "GLBMUJPSI_TEST7TEV_TABLE"
process.MuonPerformanceESProducer_GlobalMuon2.WorkingPointName = "GLBMUJPSI_TEST7TEV_WP"
process.MuonPerformanceESProducer_GlobalMuon2.ComponentName = "GlobalMuonFromTrackerTrackJpsi"
process.MuonPerformanceESProducer_TrackerTrackMuon2.PayloadName = "TRKEFFMUJPSI_TEST7TEV_TABLE"
process.MuonPerformanceESProducer_TrackerTrackMuon2.WorkingPointName = "TRKEFFMUJPSI_TEST7TEV_WP"
process.MuonPerformanceESProducer_TrackerTrackMuon2.ComponentName = "TrackerTrackFromStandaloneMuonJpsi"
process.MuonPerformanceESProducer_TriggerMuon2.PayloadName = "TRGMUJPSI_TEST7TEV_TABLE"
process.MuonPerformanceESProducer_TriggerMuon2.WorkingPointName = "TRGMUJPSI_TEST7TEV_WP"
process.MuonPerformanceESProducer_TriggerMuon2.ComponentName = "TriggerMuonFromGlobalMuonJpsi"

# PAG-specific selections
process.MuonPerformanceESProducer_TriggerMuon3.PayloadName = "TRGMUJPSI_JPSIANAL_OCTXTEST7TEV_TABLE"
process.MuonPerformanceESProducer_TriggerMuon3.WorkingPointName = "TRGMUJPSI_JPSIANAL_OCTXTEST7TEV_WP"
process.MuonPerformanceESProducer_TriggerMuon3.ComponentName = "TriggerMuonFromGlobalMuonJpsi_JpsiAnal"
process.MuonPerformanceESProducer_Muon1.PayloadName = "MUJPSI_JPSIGLBANAL_OCTXTEST7TEV_TABLE"
process.MuonPerformanceESProducer_Muon1.WorkingPointName = "MUJPSI_JPSIGLBANAL_OCTXTEST7TEV_WP"
process.MuonPerformanceESProducer_Muon1.ComponentName = "MuonFromTrackerTrackJpsi_JpsiGlbAnal"
process.MuonPerformanceESProducer_Muon2.PayloadName = "MUJPSI_JPSITKMANAL_OCTXTEST7TEV_TABLE"
process.MuonPerformanceESProducer_Muon2.WorkingPointName = "MUJPSI_JPSITKMANAL_OCTXTEST7TEV_WP"
process.MuonPerformanceESProducer_Muon2.ComponentName = "MuonFromTrackerTrackJpsi_JpsiTkMAnal"
process.MuonPerformanceESProducer_Muon3.PayloadName = "MUJPSI_JPSIPLUSMUANAL_OCTXTEST7TEV_TABLE"
process.MuonPerformanceESProducer_Muon3.WorkingPointName = "MUJPSI_JPSIPLUSMUANAL_OCTXTEST7TEV_WP"
process.MuonPerformanceESProducer_Muon3.ComponentName = "MuonFromTrackerTrackJpsi_JpsiPlusMuAnal"
process.MuonPerformanceESProducer_Muon4.PayloadName = "MUJPSI_BEXCLANAL_OCTXTEST7TEV_TABLE"
process.MuonPerformanceESProducer_Muon4.WorkingPointName = "MUJPSI_BEXCLANAL_OCTXTEST7TEV_WP"
process.MuonPerformanceESProducer_Muon4.ComponentName = "MuonFromTrackerTrackJpsi_BExclAnal"
process.MuonPerformanceESProducer_TriggerMuon4.PayloadName = "TRGMUZ_WMUNUANAL_OCTXTEST7TEV_TABLE"
process.MuonPerformanceESProducer_TriggerMuon4.WorkingPointName = "TRGMUZ_WMUNUANAL_OCTXTEST7TEV_WP"
process.MuonPerformanceESProducer_TriggerMuon4.ComponentName = "TriggerMuonFromGlobalMuonZ_WmunuAnal"
process.MuonPerformanceESProducer_Muon5.PayloadName = "MUZ_WMUNUANAL_OCTXTEST7TEV_TABLE"
process.MuonPerformanceESProducer_Muon5.WorkingPointName = "MUZ_WMUNUANAL_OCTXTEST7TEV_WP"
process.MuonPerformanceESProducer_Muon5.ComponentName = "MuonFromTrackerTrackZ_WmunuAnal"


process.selectedPatMuonsWithEff = cms.EDAnalyzer('MuTestPAT',
                                                    AlgoNames = cms.vstring(
    'GlobalMuonFromTrackerTrackJpsi'
    ))

# JH- testing PAT stuff
process.load("PhysicsTools.PatAlgos.producersLayer1.muonProducer_cff")
process.load("PhysicsTools.PatAlgos.selectionLayer1.muonSelector_cfi")
process.output = cms.OutputModule("PoolOutputModule",
                                  outputCommands = cms.untracked.vstring('drop *',
                                                                         'keep *_selectedPatMuonsWithEff_*_*'),
                                  fileName = cms.untracked.string('file:/tmp/jjhollar/jpsimm.pattuple.root'))
#                                  dataset = cms.untracked.PSet(
#    dataTier = cms.untracked.string('PAT'),
#    filterName = cms.untracked.string('')
#    )
#                                  )


process.p = cms.Path(
    #    process.patDefaultSequence +
    process.makePatMuons *
    process.selectedPatMuons +
    process.selectedPatMuonsWithEff +
    process.output)

#print process.dumpPython()

#


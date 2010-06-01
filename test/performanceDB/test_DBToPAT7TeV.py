import FWCore.ParameterSet.Config as cms

process = cms.Process("MUONEFF")

from PhysicsTools.PatAlgos.patTemplate_cfg import *

process.load("CondCore.DBCommon.CondDBCommon_cfi")

process.CondDBCommon.connect = 'sqlite_file:MuonPhysicsPerformance7TeV.db'

process.load ("MuonAnalysis.TagAndProbe.MuonPerformanceESProducer_cfi")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
        '/store/relval/CMSSW_3_6_0/RelValJpsiMM/GEN-SIM-RECO/START36_V4-v1/0013/D4B634F3-8149-DF11-9056-002618943939.root'
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
    ),
    cms.PSet(
    record = cms.string('PerformanceWPRecord'),
    tag = cms.string('TRGMUJPSI_TEST7TEV_WP'),
    label = cms.untracked.string('TRGMUJPSI_TEST7TEV_WP')
    ),                                              
    cms.PSet(
    record = cms.string('PerformancePayloadRecord'),
    tag = cms.string('TRGMUJPSI_TEST7TEV_TABLE'),
    label = cms.untracked.string('TRGMUJPSI_TEST7TEV_TABLE')
    )))


process.selectedPatMuonsWithEff = cms.EDAnalyzer('MuTestPAT',
                                    AlgoNames = cms.vstring(
                                                     'GlobalMuonFromTrackerTrackJpsi',
                                                     'TriggerMuonFromGlobalMuonJpsi'
                                    ))

process.output = cms.OutputModule("PoolOutputModule",
                                  outputCommands = cms.untracked.vstring('drop *',
                                                                         'keep *_selectedPatMuonsWithEff_*_*'),
                                  fileName = cms.untracked.string('file:/tmp/jjhollar/jpsimm.pattuple.root'))

process.p = cms.Path(
    process.makePatMuons *
    process.selectedPatMuons +
    process.selectedPatMuonsWithEff +
    process.output)

#print process.dumpPython()




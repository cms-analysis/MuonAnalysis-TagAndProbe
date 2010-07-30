import FWCore.ParameterSet.Config as cms

process = cms.Process("MUONEFF")

from PhysicsTools.PatAlgos.patTemplate_cfg import *

process.load("MuonAnalysis.TagAndProbe.MuonPerformanceESSource_cfi")

process.poolDBESSource.connect = 'sqlite_file:MuonPhysicsPerformance7TeV.db'

process.load ("MuonAnalysis.TagAndProbe.MuonPerformanceESProducer_cfi")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
        '/store/relval/CMSSW_3_6_0/RelValJpsiMM/GEN-SIM-RECO/START36_V4-v1/0013/D4B634F3-8149-DF11-9056-002618943939.root'
    )
                            )

process.selectedPatMuonsWithEff = cms.EDAnalyzer('MuTestPAT',
                                    UseAbsEtaVals = cms.bool(True),             
                                    AlgoNames = cms.vstring(
                                            'GlobalMuon_Data_CaloMuonProbe_JPsi',
                                            'HLT_Mu3_Data_CaloMuonProbe_JPsi'
                                    ))

process.out.fileName = 'file:jpsimm.pattuple.root'
process.out.outputCommands = ['drop *','keep *_selectedPatMuonsWithEff_*_*']

process.p = cms.Path(
    process.makePatMuons *
    process.selectedPatMuons +
    process.selectedPatMuonsWithEff)

#print process.dumpPython()




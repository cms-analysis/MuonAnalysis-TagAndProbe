import FWCore.ParameterSet.Config as cms

process = cms.Process("ANALYSIS")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
    'file:jpsimm.pattuple.root'
    )
                            )

process.demo2 = cms.EDAnalyzer('PATTupleReadTest')

process.p = cms.Path(
    process.demo2)


import FWCore.ParameterSet.Config as cms

process = cms.Process("myprocess")
process.load("CondCore.DBCommon.CondDBCommon_cfi")

process.CondDBCommon.connect = 'sqlite_file:MuonPhysicsPerformance.db'


process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1)
)
process.source = cms.Source("EmptySource")

process.PoolDBOutputService = cms.Service("PoolDBOutputService",
    process.CondDBCommon,
    timetype = cms.untracked.string("runnumber"),                                          
    toPut = cms.VPSet(cms.PSet(
        record = cms.string('SAMU_TABLE'),
        tag = cms.string('SAMU_TABLE'),
        label = cms.string('SAMU_TABLE')
    ),
cms.PSet(
        record = cms.string('SAMU_WP'),
        tag = cms.string('SAMU_WP'),
        label = cms.string('SAMU_WP')
    ))
                      
)

#process.mywriter = cms.EDFilter("PhysicsPerformanceDBWriterFromFile_WPandPayload",
process.mywriter = cms.EDFilter("PhysicsPerformanceDBWriterFromFile_WPandPayload_IOV",
                                inputTxtFile = cms.untracked.string('exampleStandaloneMuon.txt'),
                                RecordPayload = cms.untracked.string('SAMU_TABLE'),
                                RecordWP = cms.untracked.string('SAMU_WP'),
                                IOVBegin = cms.uint64(1),
                                IOVEnd = cms.uint64(50)
                                
#                                IOVBegin = cms.uint64(1),
#                                IOVEnd = cms.uint64(4294967295)
                                )



process.p = cms.Path(process.mywriter)



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
    toPut = cms.VPSet(
    cms.PSet(
    record = cms.string('GLBMU_TABLE'),
    tag = cms.string('GLBMU_TABLE'),
    label = cms.string('GLBMU_TABLE')
    ),
    cms.PSet(
    record = cms.string('GLBMU_WP'),
    tag = cms.string('GLBMU_WP'),
    label = cms.string('GLBMU_WP')
    ),
    cms.PSet(
    record = cms.string('TRKEFFMU_TABLE'),
    tag = cms.string('TRKEFFMU_TABLE'),
    label = cms.string('TRKEFFMU_TABLE')
    ),
    cms.PSet(
    record = cms.string('TRKEFFMU_WP'),
    tag = cms.string('TRKEFFMU_WP'),
    label = cms.string('TRKEFFMU_WP')
    ),
    cms.PSet(
    record = cms.string('TRGMU_TABLE'),
    tag = cms.string('TRGMU_TABLE'),
    label = cms.string('TRGMU_TABLE')
    ),
    cms.PSet(
    record = cms.string('TRGMU_WP'),
    tag = cms.string('TRGMU_WP'),
    label = cms.string('TRGMU_WP')
    ),    
    cms.PSet(
    record = cms.string('GLBMUCAL_TABLE'),
    tag = cms.string('GLBMUCAL_TABLE'),
    label = cms.string('GLBMUCAL_TABLE')
    ),
    cms.PSet(
    record = cms.string('GLBMUCAL_WP'),
    tag = cms.string('GLBMUCAL_WP'),
    label = cms.string('GLBMUCAL_WP')
    )
    )
                                          )
                      



process.mywriter = cms.EDFilter("PhysicsPerformanceDBWriterFromTPHist",

                                # For each table to be loaded, set the name of the input T/P file, histogram, algorithm, and cut (if any)
                                inputHistoFiles = cms.vstring('../jpsi/fit_result_GlbFromTk.root',
                                                              '../jpsi/fit_result_TkFromSta.root',
                                                              '../jpsi/fit_result_HltFromGlb.root',
                                                              '../jpsi/fit_result_GlbFromCal.root'),                                

                                inputHistogramNames = cms.vstring('fit_eff_Pt_Eta',
                                                                  'fit_eff_Pt_Eta',
                                                                  'fit_eff_Pt_Eta',
                                                                  'fit_eff_Pt_Eta'),

                                inputAlgorithmNames = cms.vstring('GlobalMuonFromTrackerTrackJpsi',
                                                                  'TrackerTrackFromStandaloneMuonJpsi',
                                                                  'TriggerMuonFromGlobalMuonJpsi',
                                                                  'GlobalMuonFromCaloMuonJpsi'), 

                                inputDiscriminatorCuts = cms.vdouble(0.0,
                                                                     0.0,
                                                                     0.0,
                                                                     0.0),

                                # For each table to be loaded, set the payload and working point record names as
                                # defined above in the PoolDBOutputService
                                RecordPayloads = cms.vstring('GLBMU_TABLE','TRKEFFMU_TABLE','TRGMU_TABLE','GLBMUCAL_TABLE'),
                                RecordWPs = cms.vstring('GLBMU_WP','TRKEFFMU_WP','TRGMU_WP','GLBMUCAL_WP'),

                                # Set the type of data to be stored, the binning variables, and the IOV
                                # These are currently assumed to be the same for all tables loaded in a single job
                                inputResultTypes = cms.vstring('efficiency','efficiency_symerr'),
                                inputBinningVariables = cms.vstring('pt','eta'),                                
                                IOVBegin = cms.uint64(1),
                                IOVEnd = cms.uint64(50)
                                )


process.p = cms.Path(process.mywriter)



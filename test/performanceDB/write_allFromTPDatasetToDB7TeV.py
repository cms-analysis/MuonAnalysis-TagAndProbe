import FWCore.ParameterSet.Config as cms

process = cms.Process("myprocess")
process.load("CondCore.DBCommon.CondDBCommon_cfi")

process.CondDBCommon.connect = 'sqlite_file:MuonPhysicsPerformance7TeV.db'


process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1)
)
process.source = cms.Source("EmptySource")

process.PoolDBOutputService = cms.Service("PoolDBOutputService",
    process.CondDBCommon,
    timetype = cms.untracked.string("runnumber"),                                          
    toPut = cms.VPSet(
    cms.PSet(
    record = cms.string('GLBMUJPSI_TEST7TEV_TABLE'),
    tag = cms.string('GLBMUJPSI_TEST7TEV_TABLE'),
    label = cms.string('GLBMUJPSI_TEST7TEV_TABLE')
    ),
    cms.PSet(
    record = cms.string('GLBMUJPSI_TEST7TEV_WP'),
    tag = cms.string('GLBMUJPSI_TEST7TEV_WP'),
    label = cms.string('GLBMUJPSI_TEST7TEV_WP')
    ),
    cms.PSet(
    record = cms.string('TRKEFFMUJPSI_TEST7TEV_TABLE'),
    tag = cms.string('TRKEFFMUJPSI_TEST7TEV_TABLE'),
    label = cms.string('TRKEFFMUJPSI_TEST7TEV_TABLE')
    ),
    cms.PSet(
    record = cms.string('TRKEFFMUJPSI_TEST7TEV_WP'),
    tag = cms.string('TRKEFFMUJPSI_TEST7TEV_WP'),
    label = cms.string('TRKEFFMUJPSI_TEST7TEV_WP')
    ),
    cms.PSet(
    record = cms.string('TRGMUJPSI_TEST7TEV_TABLE'),
    tag = cms.string('TRGMUJPSI_TEST7TEV_TABLE'),
    label = cms.string('TRGMUJPSI_TEST7TEV_TABLE')
    ),
    cms.PSet(
    record = cms.string('TRGMUJPSI_TEST7TEV_WP'),
    tag = cms.string('TRGMUJPSI_TEST7TEV_WP'),
    label = cms.string('TRGMUJPSI_TEST7TEV_WP')
    ),    
    cms.PSet(
    record = cms.string('GLBMUJPSICAL_TEST7TEV_TABLE'),
    tag = cms.string('GLBMUJPSICAL_TEST7TEV_TABLE'),
    label = cms.string('GLBMUJPSICAL_TEST7TEV_TABLE')
    ),
    cms.PSet(
    record = cms.string('GLBMUJPSICAL_TEST7TEV_WP'),
    tag = cms.string('GLBMUJPSICAL_TEST7TEV_WP'),
    label = cms.string('GLBMUJPSICAL_TEST7TEV_WP')
    ),
    cms.PSet(
    record = cms.string('GLBMUZ_TEST7TEV_TABLE'),
    tag = cms.string('GLBMUZ_TEST7TEV_TABLE'),
    label = cms.string('GLBMUZ_TEST7TEV_TABLE')
    ),
    cms.PSet(
    record = cms.string('GLBMUZ_TEST7TEV_WP'),
    tag = cms.string('GLBMUZ_TEST7TEV_WP'),
    label = cms.string('GLBMUZ_TEST7TEV_WP')
    ),
    cms.PSet(
    record = cms.string('TRKEFFMUZ_TEST7TEV_TABLE'),
    tag = cms.string('TRKEFFMUZ_TEST7TEV_TABLE'),
    label = cms.string('TRKEFFMUZ_TEST7TEV_TABLE')
    ),
    cms.PSet(
    record = cms.string('TRKEFFMUZ_TEST7TEV_WP'),
    tag = cms.string('TRKEFFMUZ_TEST7TEV_WP'),
    label = cms.string('TRKEFFMUZ_TEST7TEV_WP')
    ),
    cms.PSet(
    record = cms.string('TRGMUZ_TEST7TEV_TABLE'),
    tag = cms.string('TRGMUZ_TEST7TEV_TABLE'),
    label = cms.string('TRGMUZ_TEST7TEV_TABLE')
    ),
    cms.PSet(
    record = cms.string('TRGMUZ_TEST7TEV_WP'),
    tag = cms.string('TRGMUZ_TEST7TEV_WP'),
    label = cms.string('TRGMUZ_TEST7TEV_WP')
    ),    
    cms.PSet(
    record = cms.string('GLBMUZCAL_TEST7TEV_TABLE'),
    tag = cms.string('GLBMUZCAL_TEST7TEV_TABLE'),
    label = cms.string('GLBMUZCAL_TEST7TEV_TABLE')
    ),
    cms.PSet(
    record = cms.string('GLBMUZCAL_TEST7TEV_WP'),
    tag = cms.string('GLBMUZCAL_TEST7TEV_WP'),
    label = cms.string('GLBMUZCAL_TEST7TEV_WP')
    )))
                      



process.mywriter = cms.EDFilter("PhysicsPerformanceDBWriterFromTPDataset",

                                # For each table to be loaded, set the name of the input T/P file, histogram, algorithm, and cut (if any)
                                inputHistoFiles = cms.vstring('TnP_MuonID_data_all.root',
                                                              'TnP_Trigger_data_all.root'),
                                
                                inputDatasetNames = cms.vstring('histoMuFromTk/POG_Glb_pt_eta/fit_eff',
                                                                'histoTrigger/HLTMu3_pt_eta/fit_eff'),
                                
                                inputAlgorithmNames = cms.vstring('GlobalMuonFromTrackerTrackJpsi',
                                                                  'TriggerMuonFromGlobalMuonJpsi'),
                                
                                inputDiscriminatorCuts = cms.vdouble(0.0,
                                                                     0.0),
                                
                                
                                # For each table to be loaded, set the payload and working point record names as
                                # defined above in the PoolDBOutputService
                                RecordPayloads = cms.vstring('GLBMUJPSI_TEST7TEV_TABLE',
                                                             'TRGMUJPSI_TEST7TEV_TABLE'),
                                
                                RecordWPs = cms.vstring('GLBMUJPSI_TEST7TEV_WP',
                                                        'TRGMUJPSI_TEST7TEV_WP'),

                                # Set the type of data to be stored, the binning variables, and the IOV
                                # These are currently assumed to be the same for all tables loaded in a single job
                                inputResultTypes = cms.vstring('efficiency','efficiency_symerr'),
                                inputBinningVariables = cms.vstring('pt','eta'),                                
                                IOVBegin = cms.uint64(1),
                                IOVEnd = cms.uint64(50)
                                )


process.p = cms.Path(process.mywriter)



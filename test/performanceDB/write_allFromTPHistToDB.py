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
    record = cms.string('GLBMUJPSI_TABLE'),
    tag = cms.string('GLBMUJPSI_TABLE'),
    label = cms.string('GLBMUJPSI_TABLE')
    ),
    cms.PSet(
    record = cms.string('GLBMUJPSI_WP'),
    tag = cms.string('GLBMUJPSI_WP'),
    label = cms.string('GLBMUJPSI_WP')
    ),
    cms.PSet(
    record = cms.string('TRKEFFMUJPSI_TABLE'),
    tag = cms.string('TRKEFFMUJPSI_TABLE'),
    label = cms.string('TRKEFFMUJPSI_TABLE')
    ),
    cms.PSet(
    record = cms.string('TRKEFFMUJPSI_WP'),
    tag = cms.string('TRKEFFMUJPSI_WP'),
    label = cms.string('TRKEFFMUJPSI_WP')
    ),
    cms.PSet(
    record = cms.string('TRGMUJPSI_TABLE'),
    tag = cms.string('TRGMUJPSI_TABLE'),
    label = cms.string('TRGMUJPSI_TABLE')
    ),
    cms.PSet(
    record = cms.string('TRGMUJPSI_WP'),
    tag = cms.string('TRGMUJPSI_WP'),
    label = cms.string('TRGMUJPSI_WP')
    ),    
    cms.PSet(
    record = cms.string('GLBMUJPSICAL_TABLE'),
    tag = cms.string('GLBMUJPSICAL_TABLE'),
    label = cms.string('GLBMUJPSICAL_TABLE')
    ),
    cms.PSet(
    record = cms.string('GLBMUJPSICAL_WP'),
    tag = cms.string('GLBMUJPSICAL_WP'),
    label = cms.string('GLBMUJPSICAL_WP')
    ),
    cms.PSet(
    record = cms.string('GLBMUZ_TABLE'),
    tag = cms.string('GLBMUZ_TABLE'),
    label = cms.string('GLBMUZ_TABLE')
    ),
    cms.PSet(
    record = cms.string('GLBMUZ_WP'),
    tag = cms.string('GLBMUZ_WP'),
    label = cms.string('GLBMUZ_WP')
    ),
    cms.PSet(
    record = cms.string('TRKEFFMUZ_TABLE'),
    tag = cms.string('TRKEFFMUZ_TABLE'),
    label = cms.string('TRKEFFMUZ_TABLE')
    ),
    cms.PSet(
    record = cms.string('TRKEFFMUZ_WP'),
    tag = cms.string('TRKEFFMUZ_WP'),
    label = cms.string('TRKEFFMUZ_WP')
    ),
    cms.PSet(
    record = cms.string('TRGMUZ_TABLE'),
    tag = cms.string('TRGMUZ_TABLE'),
    label = cms.string('TRGMUZ_TABLE')
    ),
    cms.PSet(
    record = cms.string('TRGMUZ_WP'),
    tag = cms.string('TRGMUZ_WP'),
    label = cms.string('TRGMUZ_WP')
    ),    
    cms.PSet(
    record = cms.string('GLBMUZCAL_TABLE'),
    tag = cms.string('GLBMUZCAL_TABLE'),
    label = cms.string('GLBMUZCAL_TABLE')
    ),
    cms.PSet(
    record = cms.string('GLBMUZCAL_WP'),
    tag = cms.string('GLBMUZCAL_WP'),
    label = cms.string('GLBMUZCAL_WP')
    )
    
    )
                                          )
                      



process.mywriter = cms.EDFilter("PhysicsPerformanceDBWriterFromTPHist",

                                # For each table to be loaded, set the name of the input T/P file, histogram, algorithm, and cut (if any)
                                inputHistoFiles = cms.vstring('../jpsi/fit_result_GlbFromTk.root',
                                                              '../jpsi/fit_result_TkFromSta.root',
                                                              '../jpsi/fit_result_HltFromGlb.root',
                                                              '../jpsi/fit_result_GlbFromCal.root',
                                                              '../zmumu/fit_result_GlbFromTk.root',
                                                              '../zmumu/fit_result_TkFromSta.root',
                                                              '../zmumu/fit_result_HltFromGlb.root',
                                                              '../zmumu/fit_result_GlbFromCal.root'
                                                              ),                                

                                inputHistogramNames = cms.vstring('fit_eff_Pt_Eta',
                                                                  'fit_eff_Pt_Eta',
                                                                  'fit_eff_Pt_Eta',
                                                                  'fit_eff_Pt_Eta',
                                                                  'fit_eff_Pt_Eta',
                                                                  'fit_eff_Pt_Eta',
                                                                  'fit_eff_Pt_Eta',
                                                                  'fit_eff_Pt_Eta'),
                                
                                inputAlgorithmNames = cms.vstring('GlobalMuonFromTrackerTrackJpsi',
                                                                  'TrackerTrackFromStandaloneMuonJpsi',
                                                                  'TriggerMuonFromGlobalMuonJpsi',
                                                                  'GlobalMuonFromCaloMuonJpsi',
                                                                  'GlobalMuonFromTrackerTrackZ',
                                                                  'TrackerTrackFromStandaloneMuonZ',
                                                                  'TriggerMuonFromGlobalMuonZ',
                                                                  'GlobalMuonFromCaloMuonZ'), 

                                inputDiscriminatorCuts = cms.vdouble(0.0,
                                                                     0.0,
                                                                     0.0,
                                                                     0.0,
                                                                     0.0,
                                                                     0.0,
                                                                     0.0,
                                                                     0.0),
                                
                                # For each table to be loaded, set the payload and working point record names as
                                # defined above in the PoolDBOutputService
                                RecordPayloads = cms.vstring('GLBMUJPSI_TABLE','TRKEFFMUJPSI_TABLE','TRGMUJPSI_TABLE','GLBMUJPSICAL_TABLE',
                                                             'GLBMUZ_TABLE','TRKEFFMUZ_TABLE','TRGMUZ_TABLE','GLBMUZCAL_TABLE'),
                                RecordWPs = cms.vstring('GLBMUJPSI_WP','TRKEFFMUJPSI_WP','TRGMUJPSI_WP','GLBMUJPSICAL_WP',
                                                        'GLBMUZ_WP','TRKEFFMUZ_WP','TRGMUZ_WP','GLBMUZCAL_WP'),

                                # Set the type of data to be stored, the binning variables, and the IOV
                                # These are currently assumed to be the same for all tables loaded in a single job
                                inputResultTypes = cms.vstring('efficiency','efficiency_symerr'),
                                inputBinningVariables = cms.vstring('pt','eta'),                                
                                IOVBegin = cms.uint64(1),
                                IOVEnd = cms.uint64(50)
                                )


process.p = cms.Path(process.mywriter)



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
    record = cms.string('GLBMU_DATA_CALO_JPSI_TABLE'),
    tag = cms.string('GLBMU_DATA_CALO_JPSI_TABLE'),
    label = cms.string('GLBMU_DATA_CALO_JPSI_TABLE')
    ),
    cms.PSet(
    record = cms.string('GLBMU_DATA_CALO_JPSI_WP'),
    tag = cms.string('GLBMU_DATA_CALO_JPSI_WP'),
    label = cms.string('GLBMU_DATA_CALO_JPSI_WP')
    ),
    cms.PSet(
    record = cms.string('TRK_DATA_SA_JPSI_TABLE'),
    tag = cms.string('TRK_DATA_SA_JPSI_TABLE'),
    label = cms.string('TRK_DATA_SA_JPSI_TABLE')
    ),
    cms.PSet(
    record = cms.string('TRK_DATA_SA_JPSI_WP'),
    tag = cms.string('TRK_DATA_SA_JPSI_WP'),
    label = cms.string('TRK_DATA_SA_JPSI_WP')
    ),
    cms.PSet(
    record = cms.string('HLTMU3_DATA_CALO_JPSI_TABLE'),
    tag = cms.string('HLTMU3_DATA_CALO_JPSI_TABLE'),
    label = cms.string('HLTMU3_DATA_CALO_JPSI_TABLE')
    ),
    cms.PSet(
    record = cms.string('HLTMU3_DATA_CALO_JPSI_WP'),
    tag = cms.string('HLTMU3_DATA_CALO_JPSI_WP'),
    label = cms.string('HLTMU3_DATA_CALO_JPSI_WP')
    ),    
    cms.PSet(
    record = cms.string('TRKMULSAT_DATA_CALO_JPSI_TABLE'),
    tag = cms.string('TRKMULSAT_DATA_CALO_JPSI_TABLE'),
    label = cms.string('TRKMULSAT_DATA_CALO_JPSI_TABLE'),
    ),
    cms.PSet(
    record = cms.string('TRKMULSAT_DATA_CALO_JPSI_WP'),
    tag = cms.string('TRKMULSAT_DATA_CALO_JPSI_WP'),
    label = cms.string('TRKMULSAT_DATA_CALO_JPSI_WP'),
    ),
    cms.PSet(
    record = cms.string('HLTL1DOUBLEMUOPEN_DATA_CALO_JPSI_TABLE'),
    tag = cms.string('HLTL1DOUBLEMUOPEN_DATA_CALO_JPSI_TABLE'),
    label = cms.string('HLTL1DOUBLEMUOPEN_DATA_CALO_JPSI_TABLE')
    ),
    cms.PSet(
    record = cms.string('HLTL1DOUBLEMUOPEN_DATA_CALO_JPSI_WP'),
    tag = cms.string('HLTL1DOUBLEMUOPEN_DATA_CALO_JPSI_WP'),
    label = cms.string('HLTL1DOUBLEMUOPEN_DATA_CALO_JPSI_WP')
    ),
                                              
    cms.PSet(
    record = cms.string('LOERR_GLBMU_DATA_CALO_JPSI_TABLE'),
    tag = cms.string('LOERR_GLBMU_DATA_CALO_JPSI_TABLE'),
    label = cms.string('LOERR_GLBMU_DATA_CALO_JPSI_TABLE')
    ),
    cms.PSet(
    record = cms.string('LOERR_GLBMU_DATA_CALO_JPSI_WP'),
    tag = cms.string('LOERR_GLBMU_DATA_CALO_JPSI_WP'),
    label = cms.string('LOERR_GLBMU_DATA_CALO_JPSI_WP')
    ),
    cms.PSet(
    record = cms.string('LOERR_TRK_DATA_SA_JPSI_TABLE'),
    tag = cms.string('LOERR_TRK_DATA_SA_JPSI_TABLE'),
    label = cms.string('LOERR_TRK_DATA_SA_JPSI_TABLE')
    ),
    cms.PSet(
    record = cms.string('LOERR_TRK_DATA_SA_JPSI_WP'),
    tag = cms.string('LOERR_TRK_DATA_SA_JPSI_WP'),
    label = cms.string('LOERR_TRK_DATA_SA_JPSI_WP')
    ),
    cms.PSet(
    record = cms.string('LOERR_HLTMU3_DATA_CALO_JPSI_TABLE'),
    tag = cms.string('LOERR_HLTMU3_DATA_CALO_JPSI_TABLE'),
    label = cms.string('LOERR_HLTMU3_DATA_CALO_JPSI_TABLE')
    ),
    cms.PSet(
    record = cms.string('LOERR_HLTMU3_DATA_CALO_JPSI_WP'),
    tag = cms.string('LOERR_HLTMU3_DATA_CALO_JPSI_WP'),
    label = cms.string('LOERR_HLTMU3_DATA_CALO_JPSI_WP')
    ),    
    cms.PSet(
    record = cms.string('LOERR_TRKMULSAT_DATA_CALO_JPSI_TABLE'),
    tag = cms.string('LOERR_TRKMULSAT_DATA_CALO_JPSI_TABLE'),
    label = cms.string('LOERR_TRKMULSAT_DATA_CALO_JPSI_TABLE'),
    ),
    cms.PSet(
    record = cms.string('LOERR_TRKMULSAT_DATA_CALO_JPSI_WP'),
    tag = cms.string('LOERR_TRKMULSAT_DATA_CALO_JPSI_WP'),
    label = cms.string('LOERR_TRKMULSAT_DATA_CALO_JPSI_WP'),
    ),
    cms.PSet(
    record = cms.string('LOERR_HLTL1DOUBLEMUOPEN_DATA_CALO_JPSI_TABLE'),
    tag = cms.string('LOERR_HLTL1DOUBLEMUOPEN_DATA_CALO_JPSI_TABLE'),
    label = cms.string('LOERR_HLTL1DOUBLEMUOPEN_DATA_CALO_JPSI_TABLE')
    ),
    cms.PSet(
    record = cms.string('LOERR_HLTL1DOUBLEMUOPEN_DATA_CALO_JPSI_WP'),
    tag = cms.string('LOERR_HLTL1DOUBLEMUOPEN_DATA_CALO_JPSI_WP'),
    label = cms.string('LOERR_HLTL1DOUBLEMUOPEN_DATA_CALO_JPSI_WP')
    ),
                                              
    cms.PSet(
    record = cms.string('UPERR_GLBMU_DATA_CALO_JPSI_TABLE'),
    tag = cms.string('UPERR_GLBMU_DATA_CALO_JPSI_TABLE'),
    label = cms.string('UPERR_GLBMU_DATA_CALO_JPSI_TABLE')
    ),
    cms.PSet(
    record = cms.string('UPERR_GLBMU_DATA_CALO_JPSI_WP'),
    tag = cms.string('UPERR_GLBMU_DATA_CALO_JPSI_WP'),
    label = cms.string('UPERR_GLBMU_DATA_CALO_JPSI_WP')
    ),
    cms.PSet(
    record = cms.string('UPERR_TRK_DATA_SA_JPSI_TABLE'),
    tag = cms.string('UPERR_TRK_DATA_SA_JPSI_TABLE'),
    label = cms.string('UPERR_TRK_DATA_SA_JPSI_TABLE')
    ),
    cms.PSet(
    record = cms.string('UPERR_TRK_DATA_SA_JPSI_WP'),
    tag = cms.string('UPERR_TRK_DATA_SA_JPSI_WP'),
    label = cms.string('UPERR_TRK_DATA_SA_JPSI_WP')
    ),
    cms.PSet(
    record = cms.string('UPERR_HLTMU3_DATA_CALO_JPSI_TABLE'),
    tag = cms.string('UPERR_HLTMU3_DATA_CALO_JPSI_TABLE'),
    label = cms.string('UPERR_HLTMU3_DATA_CALO_JPSI_TABLE')
    ),
    cms.PSet(
    record = cms.string('UPERR_HLTMU3_DATA_CALO_JPSI_WP'),
    tag = cms.string('UPERR_HLTMU3_DATA_CALO_JPSI_WP'),
    label = cms.string('UPERR_HLTMU3_DATA_CALO_JPSI_WP')
    ),    
    cms.PSet(
    record = cms.string('UPERR_TRKMULSAT_DATA_CALO_JPSI_TABLE'),
    tag = cms.string('UPERR_TRKMULSAT_DATA_CALO_JPSI_TABLE'),
    label = cms.string('UPERR_TRKMULSAT_DATA_CALO_JPSI_TABLE'),
    ),
    cms.PSet(
    record = cms.string('UPERR_TRKMULSAT_DATA_CALO_JPSI_WP'),
    tag = cms.string('UPERR_TRKMULSAT_DATA_CALO_JPSI_WP'),
    label = cms.string('UPERR_TRKMULSAT_DATA_CALO_JPSI_WP'),
    ),
    cms.PSet(
    record = cms.string('UPERR_HLTL1DOUBLEMUOPEN_DATA_CALO_JPSI_TABLE'),
    tag = cms.string('UPERR_HLTL1DOUBLEMUOPEN_DATA_CALO_JPSI_TABLE'),
    label = cms.string('UPERR_HLTL1DOUBLEMUOPEN_DATA_CALO_JPSI_TABLE')
    ),                                              
    cms.PSet(
    record = cms.string('UPERR_HLTL1DOUBLEMUOPEN_DATA_CALO_JPSI_WP'),
    tag = cms.string('UPERR_HLTL1DOUBLEMUOPEN_DATA_CALO_JPSI_WP'),
    label = cms.string('UPERR_HLTL1DOUBLEMUOPEN_DATA_CALO_JPSI_WP')
    ),

    # MC
    cms.PSet(
    record = cms.string('GLBMU_MC_CALO_JPSI_TABLE'),
    tag = cms.string('GLBMU_MC_CALO_JPSI_TABLE'),
    label = cms.string('GLBMU_MC_CALO_JPSI_TABLE')
    ),
    cms.PSet(
    record = cms.string('GLBMU_MC_CALO_JPSI_WP'),
    tag = cms.string('GLBMU_MC_CALO_JPSI_WP'),
    label = cms.string('GLBMU_MC_CALO_JPSI_WP')
    ),
    cms.PSet(
    record = cms.string('TRK_MC_SA_JPSI_TABLE'),
    tag = cms.string('TRK_MC_SA_JPSI_TABLE'),
    label = cms.string('TRK_MC_SA_JPSI_TABLE')
    ),
    cms.PSet(
    record = cms.string('TRK_MC_SA_JPSI_WP'),
    tag = cms.string('TRK_MC_SA_JPSI_WP'),
    label = cms.string('TRK_MC_SA_JPSI_WP')
    ),
    cms.PSet(
    record = cms.string('HLTMU3_MC_CALO_JPSI_TABLE'),
    tag = cms.string('HLTMU3_MC_CALO_JPSI_TABLE'),
    label = cms.string('HLTMU3_MC_CALO_JPSI_TABLE')
    ),
    cms.PSet(
    record = cms.string('HLTMU3_MC_CALO_JPSI_WP'),
    tag = cms.string('HLTMU3_MC_CALO_JPSI_WP'),
    label = cms.string('HLTMU3_MC_CALO_JPSI_WP')
    ),    
    cms.PSet(
    record = cms.string('TRKMULSAT_MC_CALO_JPSI_TABLE'),
    tag = cms.string('TRKMULSAT_MC_CALO_JPSI_TABLE'),
    label = cms.string('TRKMULSAT_MC_CALO_JPSI_TABLE'),
    ),
    cms.PSet(
    record = cms.string('TRKMULSAT_MC_CALO_JPSI_WP'),
    tag = cms.string('TRKMULSAT_MC_CALO_JPSI_WP'),
    label = cms.string('TRKMULSAT_MC_CALO_JPSI_WP'),
    ),
    cms.PSet(
    record = cms.string('HLTL1DOUBLEMUOPEN_MC_CALO_JPSI_TABLE'),
    tag = cms.string('HLTL1DOUBLEMUOPEN_MC_CALO_JPSI_TABLE'),
    label = cms.string('HLTL1DOUBLEMUOPEN_MC_CALO_JPSI_TABLE')
    ),
    cms.PSet(
    record = cms.string('HLTL1DOUBLEMUOPEN_MC_CALO_JPSI_WP'),
    tag = cms.string('HLTL1DOUBLEMUOPEN_MC_CALO_JPSI_WP'),
    label = cms.string('HLTL1DOUBLEMUOPEN_MC_CALO_JPSI_WP')
    ),
                                              
    cms.PSet(
    record = cms.string('LOERR_GLBMU_MC_CALO_JPSI_TABLE'),
    tag = cms.string('LOERR_GLBMU_MC_CALO_JPSI_TABLE'),
    label = cms.string('LOERR_GLBMU_MC_CALO_JPSI_TABLE')
    ),
    cms.PSet(
    record = cms.string('LOERR_GLBMU_MC_CALO_JPSI_WP'),
    tag = cms.string('LOERR_GLBMU_MC_CALO_JPSI_WP'),
    label = cms.string('LOERR_GLBMU_MC_CALO_JPSI_WP')
    ),
    cms.PSet(
    record = cms.string('LOERR_TRK_MC_SA_JPSI_TABLE'),
    tag = cms.string('LOERR_TRK_MC_SA_JPSI_TABLE'),
    label = cms.string('LOERR_TRK_MC_SA_JPSI_TABLE')
    ),
    cms.PSet(
    record = cms.string('LOERR_TRK_MC_SA_JPSI_WP'),
    tag = cms.string('LOERR_TRK_MC_SA_JPSI_WP'),
    label = cms.string('LOERR_TRK_MC_SA_JPSI_WP')
    ),
    cms.PSet(
    record = cms.string('LOERR_HLTMU3_MC_CALO_JPSI_TABLE'),
    tag = cms.string('LOERR_HLTMU3_MC_CALO_JPSI_TABLE'),
    label = cms.string('LOERR_HLTMU3_MC_CALO_JPSI_TABLE')
    ),
    cms.PSet(
    record = cms.string('LOERR_HLTMU3_MC_CALO_JPSI_WP'),
    tag = cms.string('LOERR_HLTMU3_MC_CALO_JPSI_WP'),
    label = cms.string('LOERR_HLTMU3_MC_CALO_JPSI_WP')
    ),    
    cms.PSet(
    record = cms.string('LOERR_TRKMULSAT_MC_CALO_JPSI_TABLE'),
    tag = cms.string('LOERR_TRKMULSAT_MC_CALO_JPSI_TABLE'),
    label = cms.string('LOERR_TRKMULSAT_MC_CALO_JPSI_TABLE'),
    ),
    cms.PSet(
    record = cms.string('LOERR_TRKMULSAT_MC_CALO_JPSI_WP'),
    tag = cms.string('LOERR_TRKMULSAT_MC_CALO_JPSI_WP'),
    label = cms.string('LOERR_TRKMULSAT_MC_CALO_JPSI_WP'),
    ),
    cms.PSet(
    record = cms.string('LOERR_HLTL1DOUBLEMUOPEN_MC_CALO_JPSI_TABLE'),
    tag = cms.string('LOERR_HLTL1DOUBLEMUOPEN_MC_CALO_JPSI_TABLE'),
    label = cms.string('LOERR_HLTL1DOUBLEMUOPEN_MC_CALO_JPSI_TABLE')
    ),
    cms.PSet(
    record = cms.string('LOERR_HLTL1DOUBLEMUOPEN_MC_CALO_JPSI_WP'),
    tag = cms.string('LOERR_HLTL1DOUBLEMUOPEN_MC_CALO_JPSI_WP'),
    label = cms.string('LOERR_HLTL1DOUBLEMUOPEN_MC_CALO_JPSI_WP')
    ),
                                              
    cms.PSet(
    record = cms.string('UPERR_GLBMU_MC_CALO_JPSI_TABLE'),
    tag = cms.string('UPERR_GLBMU_MC_CALO_JPSI_TABLE'),
    label = cms.string('UPERR_GLBMU_MC_CALO_JPSI_TABLE')
    ),
    cms.PSet(
    record = cms.string('UPERR_GLBMU_MC_CALO_JPSI_WP'),
    tag = cms.string('UPERR_GLBMU_MC_CALO_JPSI_WP'),
    label = cms.string('UPERR_GLBMU_MC_CALO_JPSI_WP')
    ),
    cms.PSet(
    record = cms.string('UPERR_TRK_MC_SA_JPSI_TABLE'),
    tag = cms.string('UPERR_TRK_MC_SA_JPSI_TABLE'),
    label = cms.string('UPERR_TRK_MC_SA_JPSI_TABLE')
    ),
    cms.PSet(
    record = cms.string('UPERR_TRK_MC_SA_JPSI_WP'),
    tag = cms.string('UPERR_TRK_MC_SA_JPSI_WP'),
    label = cms.string('UPERR_TRK_MC_SA_JPSI_WP')
    ),
    cms.PSet(
    record = cms.string('UPERR_HLTMU3_MC_CALO_JPSI_TABLE'),
    tag = cms.string('UPERR_HLTMU3_MC_CALO_JPSI_TABLE'),
    label = cms.string('UPERR_HLTMU3_MC_CALO_JPSI_TABLE')
    ),
    cms.PSet(
    record = cms.string('UPERR_HLTMU3_MC_CALO_JPSI_WP'),
    tag = cms.string('UPERR_HLTMU3_MC_CALO_JPSI_WP'),
    label = cms.string('UPERR_HLTMU3_MC_CALO_JPSI_WP')
    ),    
    cms.PSet(
    record = cms.string('UPERR_TRKMULSAT_MC_CALO_JPSI_TABLE'),
    tag = cms.string('UPERR_TRKMULSAT_MC_CALO_JPSI_TABLE'),
    label = cms.string('UPERR_TRKMULSAT_MC_CALO_JPSI_TABLE'),
    ),
    cms.PSet(
    record = cms.string('UPERR_TRKMULSAT_MC_CALO_JPSI_WP'),
    tag = cms.string('UPERR_TRKMULSAT_MC_CALO_JPSI_WP'),
    label = cms.string('UPERR_TRKMULSAT_MC_CALO_JPSI_WP'),
    ),
    cms.PSet(
    record = cms.string('UPERR_HLTL1DOUBLEMUOPEN_MC_CALO_JPSI_TABLE'),
    tag = cms.string('UPERR_HLTL1DOUBLEMUOPEN_MC_CALO_JPSI_TABLE'),
    label = cms.string('UPERR_HLTL1DOUBLEMUOPEN_MC_CALO_JPSI_TABLE')
    ),                                              
    cms.PSet(
    record = cms.string('UPERR_HLTL1DOUBLEMUOPEN_MC_CALO_JPSI_WP'),
    tag = cms.string('UPERR_HLTL1DOUBLEMUOPEN_MC_CALO_JPSI_WP'),
    label = cms.string('UPERR_HLTL1DOUBLEMUOPEN_MC_CALO_JPSI_WP')
    )))



process.mywriter = cms.EDFilter("PhysicsPerformanceDBWriterFromTPDataset",

                                # For each table to be loaded, set the name of the input T/P file, histogram, algorithm, and cut (if any)
                                inputHistoFiles = cms.vstring('TnP_ICHEP_MuonID_data_all.root',
                                                              'TnP_ICHEP_MuonID_data_all.root', 
                                                              'TnP_ICHEP_Trigger_data_all.root',
                                                              'TnP_ICHEP_Trigger_data_all.root',
                                                              'TnP_ICHEP_MuonID_datalike_mc.root',
                                                              'TnP_ICHEP_MuonID_datalike_mc.root', 
                                                              'TnP_ICHEP_Trigger_datalike_mc.root',
                                                              'TnP_ICHEP_Trigger_datalike_mc.root'),
                                
                                inputDatasetNames = cms.vstring('histoMuFromCal/POG_Glb_pt_abseta/fit_eff',
                                                                'histoMuFromCal/POG_Glb_pt_abseta/fit_eff',
                                                                'histoTrigger/Cal_To_Mu3_pt/fit_eff',
                                                                'histoTrigger/Cal_To_L1DoubleMuOpen_pt/fit_eff',
                                                                'histoMuFromCal/POG_Glb_pt_abseta/fit_eff',
                                                                'histoMuFromCal/POG_Glb_pt_abseta/fit_eff',
                                                                'histoTrigger/Cal_To_Mu3_pt/fit_eff',
                                                                'histoTrigger/Cal_To_L1DoubleMuOpen_pt/fit_eff'),
                                
                                inputAlgorithmNames = cms.vstring('GlobalMuon_Data_CaloMuonProbe_JPsi',
                                                                  'TrackerMuonLSAT_Data_CaloMuonProbe_JPsi',
                                                                  'HLT_Mu3_Data_CaloMuonProbe_JPsi',
                                                                  'HLT_L1DoubleMuOpen_Data_CaloMuonProbe_JPsi',
                                                                  'GlobalMuon_Data_CaloMuonProbe_JPsi',
                                                                  'TrackerMuonLSAT_Data_CaloMuonProbe_JPsi',
                                                                  'HLT_Mu3_Data_CaloMuonProbe_JPsi',
                                                                  'HLT_L1DoubleMuOpen_Data_CaloMuonProbe_JPsi'),
                                
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
                                RecordPayloads = cms.vstring('TRKMULSAT_DATA_CALO_JPSI_TABLE',
                                                             'GLBMU_DATA_CALO_JPSI_TABLE',
                                                             'HLTMU3_DATA_CALO_JPSI_TABLE',
                                                             'HLTL1DOUBLEMUOPEN_DATA_CALO_JPSI_TABLE',
                                                             'TRKMULSAT_MC_CALO_JPSI_TABLE',
                                                             'GLBMU_MC_CALO_JPSI_TABLE',
                                                             'HLTMU3_MC_CALO_JPSI_TABLE',
                                                             'HLTL1DOUBLEMUOPEN_MC_CALO_JPSI_TABLE'),
                                
                                RecordWPs = cms.vstring('TRKMULSAT_DATA_CALO_JPSI_WP',
                                                        'GLBMU_DATA_CALO_JPSI_WP',
                                                        'HLTMU3_DATA_CALO_JPSI_WP',
                                                        'HLTL1DOUBLEMUOPEN_DATA_CALO_JPSI_WP',
                                                        'TRKMULSAT_MC_CALO_JPSI_WP',
                                                        'GLBMU_MC_CALO_JPSI_WP',
                                                        'HLTMU3_MC_CALO_JPSI_WP',
                                                        'HLTL1DOUBLEMUOPEN_MC_CALO_JPSI_WP'),

                                # Set the type of data to be stored, the binning variables, and the IOV
                                # These are currently assumed to be the same for all tables loaded in a single job
                                inputResultTypes = cms.vstring('efficiency','efficiency_symerr'),
                                inputBinningVariables = cms.vstring('pt','abseta'),                                
                                IOVBegin = cms.uint64(1),
                                IOVEnd = cms.uint64(50)
                                )


process.p = cms.Path(process.mywriter)



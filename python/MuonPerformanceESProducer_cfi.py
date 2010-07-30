import FWCore.ParameterSet.Config as cms

# Data efficiency central values
MuonPerformanceESProducer_GlobalMuon1 = cms.ESProducer("MuonPerformanceESProducer",
                                                          # this is what it makes available
                                                          ComponentName = cms.string('GlobalMuon_Data_CaloMuonProbe_Z'),
                                                          # this is where it gets the payload from
                                                          PayloadName = cms.string('GLBMU_DATA_CALO_Z_TABLE'),
                                                          WorkingPointName = cms.string('GLBMU_DATA_CALO_Z_WP')
                                                          )

MuonPerformanceESProducer_TriggerMuon1 = cms.ESProducer("MuonPerformanceESProducer",
                                                          # this is what it makes available
                                                          ComponentName = cms.string('HLT_Mu3_Data_CaloMuonProbe_Z'),
                                                          # this is where it gets the payload from
                                                          PayloadName = cms.string('HLTMU3_DATA_CALO_Z_TABLE'),
                                                          WorkingPointName = cms.string('HLTMU3_DATA_CALO_Z_WP')
                                                          )

MuonPerformanceESProducer_TrackerTrackMuon1 = cms.ESProducer("MuonPerformanceESProducer",
                                                            # this is what it makes available
                                                            ComponentName = cms.string('Track_Data_StandaloneMuonProbe_Z'),
                                                            # this is where it gets the payload from
                                                            PayloadName = cms.string('TRK_DATA_SA_Z_TABLE'),
                                                            WorkingPointName = cms.string('TRK_DATA_SA_Z_WP')
                                                            )

MuonPerformanceESProducer_GlobalMuon2 = cms.ESProducer("MuonPerformanceESProducer",
                                                       # this is what it makes available
                                                       ComponentName = cms.string('GlobalMuon_Data_CaloMuonProbe_JPsi'),
                                                       # this is where it gets the payload from
                                                       PayloadName = cms.string('GLBMU_DATA_CALO_JPSI_TABLE'),
                                                       WorkingPointName = cms.string('GLBMU_DATA_CALO_JPSI_WP')
                                                       )

MuonPerformanceESProducer_TriggerMuon2 = cms.ESProducer("MuonPerformanceESProducer",
                                                        # this is what it makes available
                                                        ComponentName = cms.string('HLT_Mu3_Data_CaloMuonProbe_JPsi'),
                                                        # this is where it gets the payload from
                                                        PayloadName = cms.string('HLTMU3_DATA_CALO_JPSI_TABLE'),
                                                        WorkingPointName = cms.string('HLTMU3_DATA_CALO_JPSI_WP')
                                                        )

MuonPerformanceESProducer_TrackerTrackMuon2 = cms.ESProducer("MuonPerformanceESProducer",
                                                             # this is what it makes available
                                                             ComponentName = cms.string('Track_Data_StandaloneMuonProbe_JPsi'),
                                                             # this is where it gets the payload from
                                                             PayloadName = cms.string('TRK_DATA_SA_JPSI_TABLE'),
                                                             WorkingPointName = cms.string('TRK_DATA_SA_JPSI_WP')
                                                             )

MuonPerformanceESProducer_TrackerMuon1 = cms.ESProducer("MuonPerformanceESProducer",
                                                        # this is what it makes available
                                                        ComponentName = cms.string('TrackerMuonLSAT_Data_CaloMuonProbe_JPsi'),
                                                        # this is where it gets the payload from
                                                        PayloadName = cms.string('TRKMULSAT_DATA_CALO_JPSI_TABLE'),
                                                        WorkingPointName = cms.string('TRKMULSAT_DATA_CALO_JPSI_WP')
                                                        )

MuonPerformanceESProducer_TriggerMuon3 = cms.ESProducer("MuonPerformanceESProducer",
                                                        # this is what it makes available
                                                        ComponentName = cms.string('HLT_L1DoubleMuOpen_Data_CaloMuonProbe_JPsi'),
                                                        # this is where it gets the payload from
                                                        PayloadName = cms.string('HLTL1DOUBLEMUOPEN_DATA_CALO_JPSI_TABLE'),
                                                        WorkingPointName = cms.string('HLTL1DOUBLEMUOPEN_DATA_CALO_JPSI_WP')
                                                        )


# Data efficiency lower errors
MuonPerformanceESProducer_GlobalMuonLoErr1 = cms.ESProducer("MuonPerformanceESProducer",
                                                          # this is what it makes available
                                                          ComponentName = cms.string('GlobalMuon_Data_CaloMuonProbe_Z_LowerError'),
                                                          # this is where it gets the payload from
                                                          PayloadName = cms.string('LOERR_GLBMU_DATA_CALO_Z_TABLE'),
                                                          WorkingPointName = cms.string('LOERR_GLBMU_DATA_CALO_Z_WP')
                                                          )

MuonPerformanceESProducer_TriggerMuonLoErr1 = cms.ESProducer("MuonPerformanceESProducer",
                                                          # this is what it makes available
                                                          ComponentName = cms.string('HLT_Mu3_Data_CaloMuonProbe_Z_LowerError'),
                                                          # this is where it gets the payload from
                                                          PayloadName = cms.string('LOERR_HLTMU3_DATA_CALO_Z_TABLE'),
                                                          WorkingPointName = cms.string('LOERR_HLTMU3_DATA_CALO_Z_WP')
                                                          )

MuonPerformanceESProducer_TrackerTrackMuonLoErr1 = cms.ESProducer("MuonPerformanceESProducer",
                                                            # this is what it makes available
                                                            ComponentName = cms.string('Track_Data_StandaloneMuonProbe_Z_LowerError'),
                                                            # this is where it gets the payload from
                                                            PayloadName = cms.string('LOERR_TRK_DATA_SA_Z_TABLE'),
                                                            WorkingPointName = cms.string('LOERR_TRK_DATA_SA_Z_WP')
                                                            )

MuonPerformanceESProducer_GlobalMuonLoErr2 = cms.ESProducer("MuonPerformanceESProducer",
                                                       # this is what it makes available
                                                       ComponentName = cms.string('GlobalMuon_Data_CaloMuonProbe_JPsi_LowerError'),
                                                       # this is where it gets the payload from
                                                       PayloadName = cms.string('LOERR_GLBMU_DATA_CALO_JPSI_TABLE'),
                                                       WorkingPointName = cms.string('LOERR_GLBMU_DATA_CALO_JPSI_WP')
                                                       )

MuonPerformanceESProducer_TriggerMuonLoErr2 = cms.ESProducer("MuonPerformanceESProducer",
                                                        # this is what it makes available
                                                        ComponentName = cms.string('HLT_Mu3_Data_CaloMuonProbe_JPsi_LowerError'),
                                                        # this is where it gets the payload from
                                                        PayloadName = cms.string('LOERR_HLTMU3_DATA_CALO_JPSI_TABLE'),
                                                        WorkingPointName = cms.string('LOERR_HLTMU3_DATA_CALO_JPSI_WP')
                                                        )

MuonPerformanceESProducer_TrackerTrackMuonLoErr2 = cms.ESProducer("MuonPerformanceESProducer",
                                                             # this is what it makes available
                                                             ComponentName = cms.string('Track_Data_StandaloneMuonProbe_JPsi_LowerError'),
                                                             # this is where it gets the payload from
                                                             PayloadName = cms.string('LOERR_TRK_DATA_SA_JPSI_TABLE'),
                                                             WorkingPointName = cms.string('LOERR_TRK_DATA_SA_JPSI_WP')
                                                             )

MuonPerformanceESProducer_TrackerMuonLoErr1 = cms.ESProducer("MuonPerformanceESProducer",
                                                             # this is what it makes available
                                                             ComponentName = cms.string('TrackerMuonLSAT_Data_CaloMuonProbe_JPsi_LowerError'),
                                                             # this is where it gets the payload from
                                                             PayloadName = cms.string('LOERR_TRKMULSAT_DATA_CALO_JPSI_TABLE'),
                                                             WorkingPointName = cms.string('LOERR_TRKMULSAT_DATA_CALO_JPSI_WP')
                                                             )

MuonPerformanceESProducer_TriggerMuonLoErr3 = cms.ESProducer("MuonPerformanceESProducer",
                                                             # this is what it makes available
                                                             ComponentName = cms.string('HLT_L1DoubleMuOpen_Data_CaloMuonProbe_JPsi_LowerError'),
                                                             # this is where it gets the payload from
                                                             PayloadName = cms.string('LOERR_HLTL1DOUBLEMUOPEN_DATA_CALO_JPSI_TABLE'),
                                                             WorkingPointName = cms.string('LOERR_HLTL1DOUBLEMUOPEN_DATA_CALO_JPSI_WP')
                                                        )


# Data efficiency upper errors
MuonPerformanceESProducer_GlobalMuonUpErr1 = cms.ESProducer("MuonPerformanceESProducer",
                                                          # this is what it makes available
                                                          ComponentName = cms.string('GlobalMuon_Data_CaloMuonProbe_Z_UpperError'),
                                                          # this is where it gets the payload from
                                                          PayloadName = cms.string('UPERR_GLBMU_DATA_CALO_Z_TABLE'),
                                                          WorkingPointName = cms.string('UPERR_GLBMU_DATA_CALO_Z_WP')
                                                          )

MuonPerformanceESProducer_TriggerMuonUpErr1 = cms.ESProducer("MuonPerformanceESProducer",
                                                          # this is what it makes available
                                                          ComponentName = cms.string('HLT_Mu3_Data_CaloMuonProbe_Z_UpperError'),
                                                          # this is where it gets the payload from
                                                          PayloadName = cms.string('UPERR_HLTMU3_DATA_CALO_Z_TABLE'),
                                                          WorkingPointName = cms.string('UPERR_HLTMU3_DATA_CALO_Z_WP')
                                                          )

MuonPerformanceESProducer_TrackerTrackMuonUpErr1 = cms.ESProducer("MuonPerformanceESProducer",
                                                            # this is what it makes available
                                                            ComponentName = cms.string('Track_Data_StandaloneMuonProbe_Z_UpperError'),
                                                            # this is where it gets the payload from
                                                            PayloadName = cms.string('UPERR_TRK_DATA_SA_Z_TABLE'),
                                                            WorkingPointName = cms.string('UPERR_TRK_DATA_SA_Z_WP')
                                                            )

MuonPerformanceESProducer_GlobalMuonUpErr2 = cms.ESProducer("MuonPerformanceESProducer",
                                                       # this is what it makes available
                                                       ComponentName = cms.string('GlobalMuon_Data_CaloMuonProbe_JPsi_UpperError'),
                                                       # this is where it gets the payload from
                                                       PayloadName = cms.string('UPERR_GLBMU_DATA_CALO_JPSI_TABLE'),
                                                       WorkingPointName = cms.string('UPERR_GLBMU_DATA_CALO_JPSI_WP')
                                                       )

MuonPerformanceESProducer_TriggerMuonUpErr2 = cms.ESProducer("MuonPerformanceESProducer",
                                                        # this is what it makes available
                                                        ComponentName = cms.string('HLT_Mu3_Data_CaloMuonProbe_JPsi_UpperError'),
                                                        # this is where it gets the payload from
                                                        PayloadName = cms.string('UPERR_HLTMU3_DATA_CALO_JPSI_TABLE'),
                                                        WorkingPointName = cms.string('UPERR_HLTMU3_DATA_CALO_JPSI_WP')
                                                        )

MuonPerformanceESProducer_TrackerTrackMuonUpErr2 = cms.ESProducer("MuonPerformanceESProducer",
                                                             # this is what it makes available
                                                             ComponentName = cms.string('Track_Data_StandaloneMuonProbe_JPsi_UpperError'),
                                                             # this is where it gets the payload from
                                                             PayloadName = cms.string('UPERR_TRK_DATA_SA_JPSI_TABLE'),
                                                             WorkingPointName = cms.string('UPERR_TRK_DATA_SA_JPSI_WP')
                                                             )

MuonPerformanceESProducer_TrackerMuonUpErr1 = cms.ESProducer("MuonPerformanceESProducer",
                                                             # this is what it makes available
                                                             ComponentName = cms.string('TrackerMuonLSAT_Data_CaloMuonProbe_JPsi_UpperError'),
                                                             # this is where it gets the payload from
                                                             PayloadName = cms.string('UPERR_TRKMULSAT_DATA_CALO_JPSI_TABLE'),
                                                             WorkingPointName = cms.string('UPERR_TRKMULSAT_DATA_CALO_JPSI_WP')
                                                             )

MuonPerformanceESProducer_TriggerMuonUpErr3 = cms.ESProducer("MuonPerformanceESProducer",
                                                             # this is what it makes available
                                                             ComponentName = cms.string('HLT_L1DoubleMuOpen_Data_CaloMuonProbe_JPsi_UpperError'),
                                                             # this is where it gets the payload from
                                                             PayloadName = cms.string('UPERR_HLTL1DOUBLEMUOPEN_DATA_CALO_JPSI_TABLE'),
                                                             WorkingPointName = cms.string('UPERR_HLTL1DOUBLEMUOPEN_DATA_CALO_JPSI_WP')
                                                             )

# MC efficiency central values
MuonPerformanceESProducer_MCGlobalMuon1 = cms.ESProducer("MuonPerformanceESProducer",
                                                          # this is what it makes available
                                                          ComponentName = cms.string('GlobalMuon_MC_CaloMuonProbe_Z'),
                                                          # this is where it gets the payload from
                                                          PayloadName = cms.string('GLBMU_MC_CALO_Z_TABLE'),
                                                          WorkingPointName = cms.string('GLBMU_MC_CALO_Z_WP')
                                                          )

MuonPerformanceESProducer_MCTriggerMuon1 = cms.ESProducer("MuonPerformanceESProducer",
                                                          # this is what it makes available
                                                          ComponentName = cms.string('HLT_Mu3_MC_CaloMuonProbe_Z'),
                                                          # this is where it gets the payload from
                                                          PayloadName = cms.string('HLTMU3_MC_CALO_Z_TABLE'),
                                                          WorkingPointName = cms.string('HLTMU3_MC_CALO_Z_WP')
                                                          )

MuonPerformanceESProducer_MCTrackerTrackMuon1 = cms.ESProducer("MuonPerformanceESProducer",
                                                            # this is what it makes available
                                                            ComponentName = cms.string('Track_MC_StandaloneMuonProbe_Z'),
                                                            # this is where it gets the payload from
                                                            PayloadName = cms.string('TRK_MC_SA_Z_TABLE'),
                                                            WorkingPointName = cms.string('TRK_MC_SA_Z_WP')
                                                            )

MuonPerformanceESProducer_MCGlobalMuon2 = cms.ESProducer("MuonPerformanceESProducer",
                                                       # this is what it makes available
                                                       ComponentName = cms.string('GlobalMuon_MC_CaloMuonProbe_JPsi'),
                                                       # this is where it gets the payload from
                                                       PayloadName = cms.string('GLBMU_MC_CALO_JPSI_TABLE'),
                                                       WorkingPointName = cms.string('GLBMU_MC_CALO_JPSI_WP')
                                                       )

MuonPerformanceESProducer_MCTriggerMuon2 = cms.ESProducer("MuonPerformanceESProducer",
                                                        # this is what it makes available
                                                        ComponentName = cms.string('HLT_Mu3_MC_CaloMuonProbe_JPsi'),
                                                        # this is where it gets the payload from
                                                        PayloadName = cms.string('HLTMU3_MC_CALO_JPSI_TABLE'),
                                                        WorkingPointName = cms.string('HLTMU3_MC_CALO_JPSI_WP')
                                                        )

MuonPerformanceESProducer_MCTrackerTrackMuon2 = cms.ESProducer("MuonPerformanceESProducer",
                                                             # this is what it makes available
                                                             ComponentName = cms.string('Track_MC_StandaloneMuonProbe_JPsi'),
                                                             # this is where it gets the payload from
                                                             PayloadName = cms.string('TRK_MC_SA_JPSI_TABLE'),
                                                             WorkingPointName = cms.string('TRK_MC_SA_JPSI_WP')
                                                             )

MuonPerformanceESProducer_MCTrackerMuon1 = cms.ESProducer("MuonPerformanceESProducer",
                                                        # this is what it makes available
                                                        ComponentName = cms.string('TrackerMuonLSAT_MC_CaloMuonProbe_JPsi'),
                                                        # this is where it gets the payload from
                                                        PayloadName = cms.string('TRKMULSAT_MC_CALO_JPSI_TABLE'),
                                                        WorkingPointName = cms.string('TRKMULSAT_MC_CALO_JPSI_WP')
                                                        )

MuonPerformanceESProducer_MCTriggerMuon3 = cms.ESProducer("MuonPerformanceESProducer",
                                                        # this is what it makes available
                                                        ComponentName = cms.string('HLT_L1DoubleMuOpen_MC_CaloMuonProbe_JPsi'),
                                                        # this is where it gets the payload from
                                                        PayloadName = cms.string('HLTL1DOUBLEMUOPEN_MC_CALO_JPSI_TABLE'),
                                                        WorkingPointName = cms.string('HLTL1DOUBLEMUOPEN_MC_CALO_JPSI_WP')
                                                        )


# MC efficiency lower errors
MuonPerformanceESProducer_MCGlobalMuonLoErr1 = cms.ESProducer("MuonPerformanceESProducer",
                                                          # this is what it makes available
                                                          ComponentName = cms.string('GlobalMuon_MC_CaloMuonProbe_Z_LowerError'),
                                                          # this is where it gets the payload from
                                                          PayloadName = cms.string('LOERR_GLBMU_MC_CALO_Z_TABLE'),
                                                          WorkingPointName = cms.string('LOERR_GLBMU_MC_CALO_Z_WP')
                                                          )

MuonPerformanceESProducer_MCTriggerMuonLoErr1 = cms.ESProducer("MuonPerformanceESProducer",
                                                          # this is what it makes available
                                                          ComponentName = cms.string('HLT_Mu3_MC_CaloMuonProbe_Z_LowerError'),
                                                          # this is where it gets the payload from
                                                          PayloadName = cms.string('LOERR_HLTMU3_MC_CALO_Z_TABLE'),
                                                          WorkingPointName = cms.string('LOERR_HLTMU3_MC_CALO_Z_WP')
                                                          )

MuonPerformanceESProducer_MCTrackerTrackMuonLoErr1 = cms.ESProducer("MuonPerformanceESProducer",
                                                            # this is what it makes available
                                                            ComponentName = cms.string('Track_MC_StandaloneMuonProbe_Z_LowerError'),
                                                            # this is where it gets the payload from
                                                            PayloadName = cms.string('LOERR_TRK_MC_SA_Z_TABLE'),
                                                            WorkingPointName = cms.string('LOERR_TRK_MC_SA_Z_WP')
                                                            )

MuonPerformanceESProducer_MCGlobalMuonLoErr2 = cms.ESProducer("MuonPerformanceESProducer",
                                                       # this is what it makes available
                                                       ComponentName = cms.string('GlobalMuon_MC_CaloMuonProbe_JPsi_LowerError'),
                                                       # this is where it gets the payload from
                                                       PayloadName = cms.string('LOERR_GLBMU_MC_CALO_JPSI_TABLE'),
                                                       WorkingPointName = cms.string('LOERR_GLBMU_MC_CALO_JPSI_WP')
                                                       )

MuonPerformanceESProducer_MCTriggerMuonLoErr2 = cms.ESProducer("MuonPerformanceESProducer",
                                                        # this is what it makes available
                                                        ComponentName = cms.string('HLT_Mu3_MC_CaloMuonProbe_JPsi_LowerError'),
                                                        # this is where it gets the payload from
                                                        PayloadName = cms.string('LOERR_HLTMU3_MC_CALO_JPSI_TABLE'),
                                                        WorkingPointName = cms.string('LOERR_HLTMU3_MC_CALO_JPSI_WP')
                                                        )

MuonPerformanceESProducer_MCTrackerTrackMuonLoErr2 = cms.ESProducer("MuonPerformanceESProducer",
                                                             # this is what it makes available
                                                             ComponentName = cms.string('Track_MC_StandaloneMuonProbe_JPsi_LowerError'),
                                                             # this is where it gets the payload from
                                                             PayloadName = cms.string('LOERR_TRK_MC_SA_JPSI_TABLE'),
                                                             WorkingPointName = cms.string('LOERR_TRK_MC_SA_JPSI_WP')
                                                             )

MuonPerformanceESProducer_MCTrackerMuonLoErr1 = cms.ESProducer("MuonPerformanceESProducer",
                                                             # this is what it makes available
                                                             ComponentName = cms.string('TrackerMuonLSAT_MC_CaloMuonProbe_JPsi_LowerError'),
                                                             # this is where it gets the payload from
                                                             PayloadName = cms.string('LOERR_TRKMULSAT_MC_CALO_JPSI_TABLE'),
                                                             WorkingPointName = cms.string('LOERR_TRKMULSAT_MC_CALO_JPSI_WP')
                                                             )

MuonPerformanceESProducer_MCTriggerMuonLoErr3 = cms.ESProducer("MuonPerformanceESProducer",
                                                             # this is what it makes available
                                                             ComponentName = cms.string('HLT_L1DoubleMuOpen_MC_CaloMuonProbe_JPsi_LowerError'),
                                                             # this is where it gets the payload from
                                                             PayloadName = cms.string('LOERR_HLTL1DOUBLEMUOPEN_MC_CALO_JPSI_TABLE'),
                                                             WorkingPointName = cms.string('LOERR_HLTL1DOUBLEMUOPEN_MC_CALO_JPSI_WP')
                                                        )


# MC efficiency upper errors
MuonPerformanceESProducer_MCGlobalMuonUpErr1 = cms.ESProducer("MuonPerformanceESProducer",
                                                          # this is what it makes available
                                                          ComponentName = cms.string('GlobalMuon_MC_CaloMuonProbe_Z_UpperError'),
                                                          # this is where it gets the payload from
                                                          PayloadName = cms.string('UPERR_GLBMU_MC_CALO_Z_TABLE'),
                                                          WorkingPointName = cms.string('UPERR_GLBMU_MC_CALO_Z_WP')
                                                          )

MuonPerformanceESProducer_MCTriggerMuonUpErr1 = cms.ESProducer("MuonPerformanceESProducer",
                                                          # this is what it makes available
                                                          ComponentName = cms.string('HLT_Mu3_MC_CaloMuonProbe_Z_UpperError'),
                                                          # this is where it gets the payload from
                                                          PayloadName = cms.string('UPERR_HLTMU3_MC_CALO_Z_TABLE'),
                                                          WorkingPointName = cms.string('UPERR_HLTMU3_MC_CALO_Z_WP')
                                                          )

MuonPerformanceESProducer_MCTrackerTrackMuonUpErr1 = cms.ESProducer("MuonPerformanceESProducer",
                                                            # this is what it makes available
                                                            ComponentName = cms.string('Track_MC_StandaloneMuonProbe_Z_UpperError'),
                                                            # this is where it gets the payload from
                                                            PayloadName = cms.string('UPERR_TRK_MC_SA_Z_TABLE'),
                                                            WorkingPointName = cms.string('UPERR_TRK_MC_SA_Z_WP')
                                                            )

MuonPerformanceESProducer_MCGlobalMuonUpErr2 = cms.ESProducer("MuonPerformanceESProducer",
                                                       # this is what it makes available
                                                       ComponentName = cms.string('GlobalMuon_MC_CaloMuonProbe_JPsi_UpperError'),
                                                       # this is where it gets the payload from
                                                       PayloadName = cms.string('UPERR_GLBMU_MC_CALO_JPSI_TABLE'),
                                                       WorkingPointName = cms.string('UPERR_GLBMU_MC_CALO_JPSI_WP')
                                                       )

MuonPerformanceESProducer_MCTriggerMuonUpErr2 = cms.ESProducer("MuonPerformanceESProducer",
                                                        # this is what it makes available
                                                        ComponentName = cms.string('HLT_Mu3_MC_CaloMuonProbe_JPsi_UpperError'),
                                                        # this is where it gets the payload from
                                                        PayloadName = cms.string('UPERR_HLTMU3_MC_CALO_JPSI_TABLE'),
                                                        WorkingPointName = cms.string('UPERR_HLTMU3_MC_CALO_JPSI_WP')
                                                        )

MuonPerformanceESProducer_MCTrackerTrackMuonUpErr2 = cms.ESProducer("MuonPerformanceESProducer",
                                                             # this is what it makes available
                                                             ComponentName = cms.string('Track_MC_StandaloneMuonProbe_JPsi_UpperError'),
                                                             # this is where it gets the payload from
                                                             PayloadName = cms.string('UPERR_TRK_MC_SA_JPSI_TABLE'),
                                                             WorkingPointName = cms.string('UPERR_TRK_MC_SA_JPSI_WP')
                                                             )

MuonPerformanceESProducer_MCTrackerMuonUpErr1 = cms.ESProducer("MuonPerformanceESProducer",
                                                             # this is what it makes available
                                                             ComponentName = cms.string('TrackerMuonLSAT_MC_CaloMuonProbe_JPsi_UpperError'),
                                                             # this is where it gets the payload from
                                                             PayloadName = cms.string('UPERR_TRKMULSAT_MC_CALO_JPSI_TABLE'),
                                                             WorkingPointName = cms.string('UPERR_TRKMULSAT_MC_CALO_JPSI_WP')
                                                             )

MuonPerformanceESProducer_MCTriggerMuonUpErr3 = cms.ESProducer("MuonPerformanceESProducer",
                                                             # this is what it makes available
                                                             ComponentName = cms.string('HLT_L1DoubleMuOpen_MC_CaloMuonProbe_JPsi_UpperError'),
                                                             # this is where it gets the payload from
                                                             PayloadName = cms.string('UPERR_HLTL1DOUBLEMUOPEN_MC_CALO_JPSI_TABLE'),
                                                             WorkingPointName = cms.string('UPERR_HLTL1DOUBLEMUOPEN_MC_CALO_JPSI_WP')
                                                             )
 

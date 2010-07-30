import FWCore.ParameterSet.Config as cms

# Efficiency central values
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


# Efficiency lower errors
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


# Efficiency upper errors
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

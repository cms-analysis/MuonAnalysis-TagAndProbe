import FWCore.ParameterSet.Config as cms

# Efficiency central values
MuonPerformanceESProducer_GlobalMuon1 = cms.ESProducer("MuonPerformanceESProducer",
                                                          # this is what it makes available
                                                          ComponentName = cms.string('GlobalMuonFromTrackerTrackZ'),
                                                          # this is where it gets the payload from
                                                          PayloadName = cms.string('GLBMUZ_TEST7TEV_TABLE'),
                                                          WorkingPointName = cms.string('GLBMUZ_TEST7TEV_WP')
                                                          )

MuonPerformanceESProducer_TriggerMuon1 = cms.ESProducer("MuonPerformanceESProducer",
                                                          # this is what it makes available
                                                          ComponentName = cms.string('TriggerMuonFromGlobalMuonZ'),
                                                          # this is where it gets the payload from
                                                          PayloadName = cms.string('TRGMUZ_TEST7TEV_TABLE'),
                                                          WorkingPointName = cms.string('TRGMUZ_TEST7TEV_WP')
                                                          )

MuonPerformanceESProducer_TrackerTrackMuon1 = cms.ESProducer("MuonPerformanceESProducer",
                                                            # this is what it makes available
                                                            ComponentName = cms.string('TrackerTrackFromStandaloneMuonZ'),
                                                            # this is where it gets the payload from
                                                            PayloadName = cms.string('TRKEFFMUZ_TEST7TEV_TABLE'),
                                                            WorkingPointName = cms.string('TRKEFFMUZ_TEST7TEV_WP')
                                                            )

MuonPerformanceESProducer_GlobalMuon2 = cms.ESProducer("MuonPerformanceESProducer",
                                                       # this is what it makes available
                                                       ComponentName = cms.string('GlobalMuonFromTrackerTrackJpsi'),
                                                       # this is where it gets the payload from
                                                       PayloadName = cms.string('GLBMUJPSI_TEST7TEV_TABLE'),
                                                       WorkingPointName = cms.string('GLBMUJPSI_TEST7TEV_WP')
                                                       )

MuonPerformanceESProducer_TriggerMuon2 = cms.ESProducer("MuonPerformanceESProducer",
                                                        # this is what it makes available
                                                        ComponentName = cms.string('TriggerMuonFromGlobalMuonJpsi'),
                                                        # this is where it gets the payload from
                                                        PayloadName = cms.string('TRGMUJPSI_TEST7TEV_TABLE'),
                                                        WorkingPointName = cms.string('TRGMUJPSI_TEST7TEV_WP')
                                                        )

MuonPerformanceESProducer_TrackerTrackMuon2 = cms.ESProducer("MuonPerformanceESProducer",
                                                             # this is what it makes available
                                                             ComponentName = cms.string('TrackerTrackFromStandaloneMuonJpsi'),
                                                             # this is where it gets the payload from
                                                             PayloadName = cms.string('TRKEFFMUJPSI_TEST7TEV_TABLE'),
                                                             WorkingPointName = cms.string('TRKEFFMUJPSI_TEST7TEV_WP')
                                                             )

# Efficiency lower errors
MuonPerformanceESProducer_GlobalMuonLoErr1 = cms.ESProducer("MuonPerformanceESProducer",
                                                          # this is what it makes available
                                                          ComponentName = cms.string('GlobalMuonFromTrackerTrackZ_LowerError'),
                                                          # this is where it gets the payload from
                                                          PayloadName = cms.string('LOERR_GLBMUZ_TEST7TEV_TABLE'),
                                                          WorkingPointName = cms.string('LOERR_GLBMUZ_TEST7TEV_WP')
                                                          )

MuonPerformanceESProducer_TriggerMuonLoErr1 = cms.ESProducer("MuonPerformanceESProducer",
                                                          # this is what it makes available
                                                          ComponentName = cms.string('TriggerMuonFromGlobalMuonZ_LowerError'),
                                                          # this is where it gets the payload from
                                                          PayloadName = cms.string('LOERR_TRGMUZ_TEST7TEV_TABLE'),
                                                          WorkingPointName = cms.string('LOERR_TRGMUZ_TEST7TEV_WP')
                                                          )

MuonPerformanceESProducer_TrackerTrackMuonLoErr1 = cms.ESProducer("MuonPerformanceESProducer",
                                                            # this is what it makes available
                                                            ComponentName = cms.string('TrackerTrackFromStandaloneMuonZ_LowerError'),
                                                            # this is where it gets the payload from
                                                            PayloadName = cms.string('LOERR_TRKEFFMUZ_TEST7TEV_TABLE'),
                                                            WorkingPointName = cms.string('LOERR_TRKEFFMUZ_TEST7TEV_WP')
                                                            )

MuonPerformanceESProducer_GlobalMuonLoErr2 = cms.ESProducer("MuonPerformanceESProducer",
                                                       # this is what it makes available
                                                       ComponentName = cms.string('GlobalMuonFromTrackerTrackJpsi_LowerError'),
                                                       # this is where it gets the payload from
                                                       PayloadName = cms.string('LOERR_GLBMUJPSI_TEST7TEV_TABLE'),
                                                       WorkingPointName = cms.string('LOERR_GLBMUJPSI_TEST7TEV_WP')
                                                       )

MuonPerformanceESProducer_TriggerMuonLoErr2 = cms.ESProducer("MuonPerformanceESProducer",
                                                        # this is what it makes available
                                                        ComponentName = cms.string('TriggerMuonFromGlobalMuonJpsi_LowerError'),
                                                        # this is where it gets the payload from
                                                        PayloadName = cms.string('LOERR_TRGMUJPSI_TEST7TEV_TABLE'),
                                                        WorkingPointName = cms.string('LOERR_TRGMUJPSI_TEST7TEV_WP')
                                                        )

MuonPerformanceESProducer_TrackerTrackMuonLoErr2 = cms.ESProducer("MuonPerformanceESProducer",
                                                             # this is what it makes available
                                                             ComponentName = cms.string('TrackerTrackFromStandaloneMuonJpsi_LowerError'),
                                                             # this is where it gets the payload from
                                                             PayloadName = cms.string('LOERR_TRKEFFMUJPSI_TEST7TEV_TABLE'),
                                                             WorkingPointName = cms.string('LOERR_TRKEFFMUJPSI_TEST7TEV_WP')
                                                             )

# Efficiency upper errors
MuonPerformanceESProducer_GlobalMuonUpErr1 = cms.ESProducer("MuonPerformanceESProducer",
                                                          # this is what it makes available
                                                          ComponentName = cms.string('GlobalMuonFromTrackerTrackZ_UpperError'),
                                                          # this is where it gets the payload from
                                                          PayloadName = cms.string('UPERR_GLBMUZ_TEST7TEV_TABLE'),
                                                          WorkingPointName = cms.string('UPERR_GLBMUZ_TEST7TEV_WP')
                                                          )

MuonPerformanceESProducer_TriggerMuonUpErr1 = cms.ESProducer("MuonPerformanceESProducer",
                                                          # this is what it makes available
                                                          ComponentName = cms.string('TriggerMuonFromGlobalMuonZ_UpperError'),
                                                          # this is where it gets the payload from
                                                          PayloadName = cms.string('UPERR_TRGMUZ_TEST7TEV_TABLE'),
                                                          WorkingPointName = cms.string('UPERR_TRGMUZ_TEST7TEV_WP')
                                                          )

MuonPerformanceESProducer_TrackerTrackMuonUpErr1 = cms.ESProducer("MuonPerformanceESProducer",
                                                            # this is what it makes available
                                                            ComponentName = cms.string('TrackerTrackFromStandaloneMuonZ_UpperError'),
                                                            # this is where it gets the payload from
                                                            PayloadName = cms.string('UPERR_TRKEFFMUZ_TEST7TEV_TABLE'),
                                                            WorkingPointName = cms.string('UPERR_TRKEFFMUZ_TEST7TEV_WP')
                                                            )

MuonPerformanceESProducer_GlobalMuonUpErr2 = cms.ESProducer("MuonPerformanceESProducer",
                                                       # this is what it makes available
                                                       ComponentName = cms.string('GlobalMuonFromTrackerTrackJpsi_UpperError'),
                                                       # this is where it gets the payload from
                                                       PayloadName = cms.string('UPERR_GLBMUJPSI_TEST7TEV_TABLE'),
                                                       WorkingPointName = cms.string('UPERR_GLBMUJPSI_TEST7TEV_WP')
                                                       )

MuonPerformanceESProducer_TriggerMuonUpErr2 = cms.ESProducer("MuonPerformanceESProducer",
                                                        # this is what it makes available
                                                        ComponentName = cms.string('TriggerMuonFromGlobalMuonJpsi_UpperError'),
                                                        # this is where it gets the payload from
                                                        PayloadName = cms.string('UPERR_TRGMUJPSI_TEST7TEV_TABLE'),
                                                        WorkingPointName = cms.string('UPERR_TRGMUJPSI_TEST7TEV_WP')
                                                        )

MuonPerformanceESProducer_TrackerTrackMuonUpErr2 = cms.ESProducer("MuonPerformanceESProducer",
                                                             # this is what it makes available
                                                             ComponentName = cms.string('TrackerTrackFromStandaloneMuonJpsi_UpperError'),
                                                             # this is where it gets the payload from
                                                             PayloadName = cms.string('UPERR_TRKEFFMUJPSI_TEST7TEV_TABLE'),
                                                             WorkingPointName = cms.string('UPERR_TRKEFFMUJPSI_TEST7TEV_WP')
                                                             )

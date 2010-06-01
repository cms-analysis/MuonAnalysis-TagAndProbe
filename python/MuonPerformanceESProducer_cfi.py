import FWCore.ParameterSet.Config as cms

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


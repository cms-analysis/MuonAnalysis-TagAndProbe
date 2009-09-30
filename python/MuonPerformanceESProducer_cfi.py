import FWCore.ParameterSet.Config as cms

MuonPerformanceESProducer_StandaloneMuon = cms.ESProducer("MuonPerformanceESProducer",
                                                          # this is what it makes available
                                                          ComponentName = cms.string('StandaloneMuon'),
                                                          # this is where it gets the payload from                                                
                                                          PayloadName = cms.string('StandaloneMuon'),
                                                          WorkingPointName = cms.string('StandaloneMuon')
                                                          )

MuonPerformanceESProducer_TrackerTrackMuon = cms.ESProducer("MuonPerformanceESProducer",
                                                            # this is what it makes available
                                                            ComponentName = cms.string('TrackerTrackMuon'),
                                                            # this is where it gets the payload from
                                                            PayloadName = cms.string('TrackerTrackMuon'),
                                                            WorkingPointName = cms.string('TrackerTrackMuon')
                                                            )


MuonPerformanceESProducer_TriggerMuon = cms.ESProducer("MuonPerformanceESProducer",
                                                       # this is what it makes available
                                                       ComponentName = cms.string('TriggerMuon'),
                                                       # this is where it gets the payload from
                                                       PayloadName = cms.string('TriggerMuon'),
                                                       WorkingPointName = cms.string('TriggerMuon')
                                                       )

import FWCore.ParameterSet.Config as cms

from CondCore.DBCommon.CondDBCommon_cfi import *
poolDBESSource = cms.ESSource("PoolDBESSource",
                              CondDBCommon,
                              toGet = cms.VPSet(
    cms.PSet(
    record = cms.string('PerformanceWPRecord'),
    tag = cms.string('GLBMU_DATA_CALO_JPSI_WP'),
    label = cms.untracked.string('GLBMU_DATA_CALO_JPSI_WP')
    ),
    cms.PSet(
    record = cms.string('PerformancePayloadRecord'),
    tag = cms.string('GLBMU_DATA_CALO_JPSI_TABLE'),
    label = cms.untracked.string('GLBMU_DATA_CALO_JPSI_TABLE')
    ),
    cms.PSet(
    record = cms.string('PerformanceWPRecord'),
    tag = cms.string('HLTMU3_DATA_CALO_JPSI_WP'),
    label = cms.untracked.string('HLTMU3_DATA_CALO_JPSI_WP')
    ),
    cms.PSet(
    record = cms.string('PerformancePayloadRecord'),
    tag = cms.string('HLTMU3_DATA_CALO_JPSI_TABLE'),
    label = cms.untracked.string('HLTMU3_DATA_CALO_JPSI_TABLE')
    ),
    cms.PSet(
    record = cms.string('PerformancePayloadRecord'),
    tag = cms.string('TRKMULSAT_DATA_CALO_JPSI_WP'),
    label = cms.untracked.string('TRKMULSAT_DATA_CALO_JPSI_WP')
    ),
    cms.PSet(
    record = cms.string('PerformancePayloadRecord'),
    tag = cms.string('TRKMULSAT_DATA_CALO_JPSI_TABLE'),
    label = cms.untracked.string('TRKMULSAT_DATA_CALO_JPSI_TABLE')
    ),
    cms.PSet(
    record = cms.string('PerformancePayloadRecord'),
    tag = cms.string('HLTL1DOUBLEMUOPEN_DATA_CALO_JPSI_WP'),
    label = cms.untracked.string('HLTL1DOUBLEMUOPEN_DATA_CALO_JPSI_WP')
    ),
    cms.PSet(
    record = cms.string('PerformancePayloadRecord'),
    tag = cms.string('HLTL1DOUBLEMUOPEN_DATA_CALO_JPSI_TABLE'),
    label = cms.untracked.string('HLTL1DOUBLEMUOPEN_DATA_CALO_JPSI_TABLE')
    ),
                                  
    cms.PSet(
    record = cms.string('PerformanceWPRecord'),
    tag = cms.string('LOERR_GLBMU_DATA_CALO_JPSI_WP'),
    label = cms.untracked.string('LOERR_GLBMU_DATA_CALO_JPSI_WP')
    ),
    cms.PSet(
    record = cms.string('PerformancePayloadRecord'),
    tag = cms.string('LOERR_GLBMU_DATA_CALO_JPSI_TABLE'),
    label = cms.untracked.string('LOERR_GLBMU_DATA_CALO_JPSI_TABLE')
    ),
    cms.PSet(
    record = cms.string('PerformanceWPRecord'),
    tag = cms.string('LOERR_HLTMU3_DATA_CALO_JPSI_WP'),
    label = cms.untracked.string('LOERR_HLTMU3_DATA_CALO_JPSI_WP')
    ),
    cms.PSet(
    record = cms.string('PerformancePayloadRecord'),
    tag = cms.string('LOERR_HLTMU3_DATA_CALO_JPSI_TABLE'),
    label = cms.untracked.string('LOERR_HLTMU3_DATA_CALO_JPSI_TABLE')
    ),
    cms.PSet(
    record = cms.string('PerformancePayloadRecord'),
    tag = cms.string('LOERR_TRKMULSAT_DATA_CALO_JPSI_WP'),
    label = cms.untracked.string('LOERR_TRKMULSAT_DATA_CALO_JPSI_WP')
    ),
    cms.PSet(
    record = cms.string('PerformancePayloadRecord'),
    tag = cms.string('LOERR_TRKMULSAT_DATA_CALO_JPSI_TABLE'),
    label = cms.untracked.string('LOERR_TRKMULSAT_DATA_CALO_JPSI_TABLE')
    ),
    cms.PSet(
    record = cms.string('PerformancePayloadRecord'),
    tag = cms.string('LOERR_HLTL1DOUBLEMUOPEN_DATA_CALO_JPSI_WP'),
    label = cms.untracked.string('LOERR_HLTL1DOUBLEMUOPEN_DATA_CALO_JPSI_WP')
    ),
    cms.PSet(
    record = cms.string('PerformancePayloadRecord'),
    tag = cms.string('LOERR_HLTL1DOUBLEMUOPEN_DATA_CALO_JPSI_TABLE'),
    label = cms.untracked.string('LOERR_HLTL1DOUBLEMUOPEN_DATA_CALO_JPSI_TABLE')
    ),

    cms.PSet(
    record = cms.string('PerformanceWPRecord'),
    tag = cms.string('UPERR_GLBMU_DATA_CALO_JPSI_WP'),
    label = cms.untracked.string('UPERR_GLBMU_DATA_CALO_JPSI_WP')
    ),
    cms.PSet(
    record = cms.string('PerformancePayloadRecord'),
    tag = cms.string('UPERR_GLBMU_DATA_CALO_JPSI_TABLE'),
    label = cms.untracked.string('UPERR_GLBMU_DATA_CALO_JPSI_TABLE')
    ),
    cms.PSet(
    record = cms.string('PerformanceWPRecord'),
    tag = cms.string('UPERR_HLTMU3_DATA_CALO_JPSI_WP'),
    label = cms.untracked.string('UPERR_HLTMU3_DATA_CALO_JPSI_WP')
    ),
    cms.PSet(
    record = cms.string('PerformancePayloadRecord'),
    tag = cms.string('UPERR_HLTMU3_DATA_CALO_JPSI_TABLE'),
    label = cms.untracked.string('UPERR_HLTMU3_DATA_CALO_JPSI_TABLE')
    ),
    cms.PSet(
    record = cms.string('PerformancePayloadRecord'),
    tag = cms.string('UPERR_TRKMULSAT_DATA_CALO_JPSI_WP'),
    label = cms.untracked.string('UPERR_TRKMULSAT_DATA_CALO_JPSI_WP')
    ),
    cms.PSet(
    record = cms.string('PerformancePayloadRecord'),
    tag = cms.string('UPERR_TRKMULSAT_DATA_CALO_JPSI_TABLE'),
    label = cms.untracked.string('UPERR_TRKMULSAT_DATA_CALO_JPSI_TABLE')
    ),
    cms.PSet(
    record = cms.string('PerformancePayloadRecord'),
    tag = cms.string('UPERR_HLTL1DOUBLEMUOPEN_DATA_CALO_JPSI_WP'),
    label = cms.untracked.string('UPERR_HLTL1DOUBLEMUOPEN_DATA_CALO_JPSI_WP')
    ),                                  
    cms.PSet(
    record = cms.string('PerformancePayloadRecord'),
    tag = cms.string('UPERR_HLTL1DOUBLEMUOPEN_DATA_CALO_JPSI_TABLE'),
    label = cms.untracked.string('UPERR_HLTL1DOUBLEMUOPEN_DATA_CALO_JPSI_TABLE')
    ),

    # MC records
    cms.PSet(
    record = cms.string('PerformanceWPRecord'),
    tag = cms.string('GLBMU_MC_CALO_JPSI_WP'),
    label = cms.untracked.string('GLBMU_MC_CALO_JPSI_WP')
    ),
    cms.PSet(
    record = cms.string('PerformancePayloadRecord'),
    tag = cms.string('GLBMU_MC_CALO_JPSI_TABLE'),
    label = cms.untracked.string('GLBMU_MC_CALO_JPSI_TABLE')
    ),
    cms.PSet(
    record = cms.string('PerformanceWPRecord'),
    tag = cms.string('HLTMU3_MC_CALO_JPSI_WP'),
    label = cms.untracked.string('HLTMU3_MC_CALO_JPSI_WP')
    ),
    cms.PSet(
    record = cms.string('PerformancePayloadRecord'),
    tag = cms.string('HLTMU3_MC_CALO_JPSI_TABLE'),
    label = cms.untracked.string('HLTMU3_MC_CALO_JPSI_TABLE')
    ),
    cms.PSet(
    record = cms.string('PerformancePayloadRecord'),
    tag = cms.string('TRKMULSAT_MC_CALO_JPSI_WP'),
    label = cms.untracked.string('TRKMULSAT_MC_CALO_JPSI_WP')
    ),
    cms.PSet(
    record = cms.string('PerformancePayloadRecord'),
    tag = cms.string('TRKMULSAT_MC_CALO_JPSI_TABLE'),
    label = cms.untracked.string('TRKMULSAT_MC_CALO_JPSI_TABLE')
    ),
    cms.PSet(
    record = cms.string('PerformancePayloadRecord'),
    tag = cms.string('HLTL1DOUBLEMUOPEN_MC_CALO_JPSI_WP'),
    label = cms.untracked.string('HLTL1DOUBLEMUOPEN_MC_CALO_JPSI_WP')
    ),
    cms.PSet(
    record = cms.string('PerformancePayloadRecord'),
    tag = cms.string('HLTL1DOUBLEMUOPEN_MC_CALO_JPSI_TABLE'),
    label = cms.untracked.string('HLTL1DOUBLEMUOPEN_MC_CALO_JPSI_TABLE')
    ),
                                  
    cms.PSet(
    record = cms.string('PerformanceWPRecord'),
    tag = cms.string('LOERR_GLBMU_MC_CALO_JPSI_WP'),
    label = cms.untracked.string('LOERR_GLBMU_MC_CALO_JPSI_WP')
    ),
    cms.PSet(
    record = cms.string('PerformancePayloadRecord'),
    tag = cms.string('LOERR_GLBMU_MC_CALO_JPSI_TABLE'),
    label = cms.untracked.string('LOERR_GLBMU_MC_CALO_JPSI_TABLE')
    ),
    cms.PSet(
    record = cms.string('PerformanceWPRecord'),
    tag = cms.string('LOERR_HLTMU3_MC_CALO_JPSI_WP'),
    label = cms.untracked.string('LOERR_HLTMU3_MC_CALO_JPSI_WP')
    ),
    cms.PSet(
    record = cms.string('PerformancePayloadRecord'),
    tag = cms.string('LOERR_HLTMU3_MC_CALO_JPSI_TABLE'),
    label = cms.untracked.string('LOERR_HLTMU3_MC_CALO_JPSI_TABLE')
    ),
    cms.PSet(
    record = cms.string('PerformancePayloadRecord'),
    tag = cms.string('LOERR_TRKMULSAT_MC_CALO_JPSI_WP'),
    label = cms.untracked.string('LOERR_TRKMULSAT_MC_CALO_JPSI_WP')
    ),
    cms.PSet(
    record = cms.string('PerformancePayloadRecord'),
    tag = cms.string('LOERR_TRKMULSAT_MC_CALO_JPSI_TABLE'),
    label = cms.untracked.string('LOERR_TRKMULSAT_MC_CALO_JPSI_TABLE')
    ),
    cms.PSet(
    record = cms.string('PerformancePayloadRecord'),
    tag = cms.string('LOERR_HLTL1DOUBLEMUOPEN_MC_CALO_JPSI_WP'),
    label = cms.untracked.string('LOERR_HLTL1DOUBLEMUOPEN_MC_CALO_JPSI_WP')
    ),
    cms.PSet(
    record = cms.string('PerformancePayloadRecord'),
    tag = cms.string('LOERR_HLTL1DOUBLEMUOPEN_MC_CALO_JPSI_TABLE'),
    label = cms.untracked.string('LOERR_HLTL1DOUBLEMUOPEN_MC_CALO_JPSI_TABLE')
    ),

    cms.PSet(
    record = cms.string('PerformanceWPRecord'),
    tag = cms.string('UPERR_GLBMU_MC_CALO_JPSI_WP'),
    label = cms.untracked.string('UPERR_GLBMU_MC_CALO_JPSI_WP')
    ),
    cms.PSet(
    record = cms.string('PerformancePayloadRecord'),
    tag = cms.string('UPERR_GLBMU_MC_CALO_JPSI_TABLE'),
    label = cms.untracked.string('UPERR_GLBMU_MC_CALO_JPSI_TABLE')
    ),
    cms.PSet(
    record = cms.string('PerformanceWPRecord'),
    tag = cms.string('UPERR_HLTMU3_MC_CALO_JPSI_WP'),
    label = cms.untracked.string('UPERR_HLTMU3_MC_CALO_JPSI_WP')
    ),
    cms.PSet(
    record = cms.string('PerformancePayloadRecord'),
    tag = cms.string('UPERR_HLTMU3_MC_CALO_JPSI_TABLE'),
    label = cms.untracked.string('UPERR_HLTMU3_MC_CALO_JPSI_TABLE')
    ),
    cms.PSet(
    record = cms.string('PerformancePayloadRecord'),
    tag = cms.string('UPERR_TRKMULSAT_MC_CALO_JPSI_WP'),
    label = cms.untracked.string('UPERR_TRKMULSAT_MC_CALO_JPSI_WP')
    ),
    cms.PSet(
    record = cms.string('PerformancePayloadRecord'),
    tag = cms.string('UPERR_TRKMULSAT_MC_CALO_JPSI_TABLE'),
    label = cms.untracked.string('UPERR_TRKMULSAT_MC_CALO_JPSI_TABLE')
    ),
    cms.PSet(
    record = cms.string('PerformancePayloadRecord'),
    tag = cms.string('UPERR_HLTL1DOUBLEMUOPEN_MC_CALO_JPSI_WP'),
    label = cms.untracked.string('UPERR_HLTL1DOUBLEMUOPEN_MC_CALO_JPSI_WP')
    ),                                  
    cms.PSet(
    record = cms.string('PerformancePayloadRecord'),
    tag = cms.string('UPERR_HLTL1DOUBLEMUOPEN_MC_CALO_JPSI_TABLE'),
    label = cms.untracked.string('UPERR_HLTL1DOUBLEMUOPEN_MC_CALO_JPSI_TABLE')                                  
    )))



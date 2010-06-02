import FWCore.ParameterSet.Config as cms

from CondCore.DBCommon.CondDBCommon_cfi import *
poolDBESSource = cms.ESSource("PoolDBESSource",
                              CondDBCommon,
                              toGet = cms.VPSet(
    cms.PSet(
    record = cms.string('PerformanceWPRecord'),
    tag = cms.string('GLBMUJPSI_TEST7TEV_WP'),
    label = cms.untracked.string('GLBMUJPSI_TEST7TEV_WP')
    ),
    cms.PSet(
    record = cms.string('PerformancePayloadRecord'),
    tag = cms.string('GLBMUJPSI_TEST7TEV_TABLE'),
    label = cms.untracked.string('GLBMUJPSI_TEST7TEV_TABLE')
    ),
    cms.PSet(
    record = cms.string('PerformanceWPRecord'),
    tag = cms.string('TRGMUJPSI_TEST7TEV_WP'),
    label = cms.untracked.string('TRGMUJPSI_TEST7TEV_WP')
    ),
    cms.PSet(
    record = cms.string('PerformancePayloadRecord'),
    tag = cms.string('TRGMUJPSI_TEST7TEV_TABLE'),
    label = cms.untracked.string('TRGMUJPSI_TEST7TEV_TABLE')
    ),
    cms.PSet(
    record = cms.string('PerformanceWPRecord'),
    tag = cms.string('LOERR_GLBMUJPSI_TEST7TEV_WP'),
    label = cms.untracked.string('LOERR_GLBMUJPSI_TEST7TEV_WP')
    ),
    cms.PSet(
    record = cms.string('PerformancePayloadRecord'),
    tag = cms.string('LOERR_GLBMUJPSI_TEST7TEV_TABLE'),
    label = cms.untracked.string('LOERR_GLBMUJPSI_TEST7TEV_TABLE')
    ),
    cms.PSet(
    record = cms.string('PerformanceWPRecord'),
    tag = cms.string('LOERR_TRGMUJPSI_TEST7TEV_WP'),
    label = cms.untracked.string('LOERR_TRGMUJPSI_TEST7TEV_WP')
    ),
    cms.PSet(
    record = cms.string('PerformancePayloadRecord'),
    tag = cms.string('LOERR_TRGMUJPSI_TEST7TEV_TABLE'),
    label = cms.untracked.string('LOERR_TRGMUJPSI_TEST7TEV_TABLE')
    ),
    cms.PSet(
    record = cms.string('PerformanceWPRecord'),
    tag = cms.string('UPERR_GLBMUJPSI_TEST7TEV_WP'),
    label = cms.untracked.string('UPERR_GLBMUJPSI_TEST7TEV_WP')
    ),
    cms.PSet(
    record = cms.string('PerformancePayloadRecord'),
    tag = cms.string('UPERR_GLBMUJPSI_TEST7TEV_TABLE'),
    label = cms.untracked.string('UPERR_GLBMUJPSI_TEST7TEV_TABLE')
    ),
    cms.PSet(
    record = cms.string('PerformanceWPRecord'),
    tag = cms.string('UPERR_TRGMUJPSI_TEST7TEV_WP'),
    label = cms.untracked.string('UPERR_TRGMUJPSI_TEST7TEV_WP')
    ),
    cms.PSet(
    record = cms.string('PerformancePayloadRecord'),
    tag = cms.string('UPERR_TRGMUJPSI_TEST7TEV_TABLE'),
    label = cms.untracked.string('UPERR_TRGMUJPSI_TEST7TEV_TABLE')
    )))



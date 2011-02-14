import FWCore.ParameterSet.Config as cms

expectedHitsComputer = cms.EDProducer("ExpectedHitsComputer",
inputColl       = cms.InputTag("muons"), 
useGsfTrack     = cms.bool(False),
objectSelection = cms.string(""),
propagator         = cms.string('PropagatorWithMaterialOpposite'),
navigationSchool   = cms.string('SimpleNavigationSchool'),
measurementTracker = cms.string(''),
)

import FWCore.ParameterSet.Config as cms

nearbyMuonsInfo = cms.EDProducer("NearbyMuonsInfo",
    # put here a collection of Composite Candidates with ShallowClones to muons
    src = cms.InputTag("tpPairs"), 
    # Configuration for the extrapolation at the muon system 
    propM1 = cms.PSet(
        useStation2 = cms.bool(False), 
        useTrack = cms.string("tracker"),
        useState = cms.string("atVertex"),  # in AOD
        useSimpleGeometry = cms.bool(True), # use just one cylinder and two planes, not all the fancy chambers
    ),
    propM2 = cms.PSet(
        useStation2 = cms.bool(True), 
        useTrack = cms.string("tracker"),
        useState = cms.string("atVertex"),  # in AOD
        useSimpleGeometry = cms.bool(True), # use just one cylinder and two planes, not all the fancy chambers
    )
)

import FWCore.ParameterSet.Config as cms

tkClusterInfo = cms.EDProducer("ComputeTkClusterInfos",
    probes = cms.InputTag("tagMuons"),
    stripClusters = cms.InputTag("siStripClusters"),
    pixelClusters = cms.InputTag("siPixelClusters"),
)
TkClusterInfoVars = cms.PSet(
    nClustersStrip = cms.InputTag("tkClusterInfo","siStripClusterCount"),
    nClustersPixel = cms.InputTag("tkClusterInfo","siPixelClusterCount"),
)

tkLogErrors = cms.EDProducer("ComputeLogErrorTotals",
    probes    = cms.InputTag("tagMuons"),
    logErrors = cms.InputTag("logErrorHarvester"),
    counters = cms.PSet(
        firstStep = cms.vstring(
            "initialStepSeeds","initialStepTracks",
        ),
        pixelSteps = cms.vstring(
            "initialStepSeeds","initialStepTracks",
            "lowPtTripletStepSeeds","lowPtTripletStepTracks",
            "pixelPairStepSeeds","pixelPairStepTracks",
            "detatchedTripletStepSeeds","detatchedTripletStepTracks",
        ), 
        anyStep = cms.vstring(
            "initialStepSeeds","initialStepTracks",
            "lowPtTripletStepSeeds","lowPtTripletStepTracks",
            "pixelPairStepSeeds","pixelPairStepTracks",
            "detatchedTripletStepSeeds","detatchedTripletStepTracks",
            "mixedTripletStepSeeds","mixedTripletStepTracks",
            "pixelLessStepSeeds","pixelLessStepTracks",
            "tobTecStepSeeds","tobTecStepTracks",
        ),
    ),
)

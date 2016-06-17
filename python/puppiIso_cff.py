import FWCore.ParameterSet.Config as cms

PuppiIsolationVariables = cms.PSet(
        muPFIsoValueCHR04PUPPI = cms.InputTag("loadPUPPIisoInValueMaps","muPFIsoValueCHR04PUPPI"),
        muPFIsoValueNHR04PUPPI = cms.InputTag("loadPUPPIisoInValueMaps","muPFIsoValueNHR04PUPPI"),
        muPFIsoValuePhR04PUPPI = cms.InputTag("loadPUPPIisoInValueMaps","muPFIsoValuePhR04PUPPI"),
        muPFIsoValueCHR04PUPPINoLep = cms.InputTag("loadPUPPIisoInValueMaps","muPFIsoValueCHR04PUPPINoLep"),
        muPFIsoValueNHR04PUPPINoLep = cms.InputTag("loadPUPPIisoInValueMaps","muPFIsoValueNHR04PUPPINoLep"),
        muPFIsoValuePhR04PUPPINoLep = cms.InputTag("loadPUPPIisoInValueMaps","muPFIsoValuePhR04PUPPINoLep"),
)

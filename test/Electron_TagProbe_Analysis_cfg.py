import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")
# keep the logging output to a nice level
process.load("FWCore.MessageService.MessageLogger_cfi")

process.load("MuonAnalysis.TagAndProbe.Electron_TagProbeEDMAnalysis_cfi")

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:output/CSA08_electron_tp_edm_ntuple_reReco.root')
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)
process.p = cms.Path(process.demo)
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.demo.TagProbeType = 1
process.demo.FitFileName = 'electron_eff_CSA08_reReco_GsfToIso_EtaPhi.root'
process.demo.NameVar1 = 'phi'
process.demo.NumBinsVar1 = 28
process.demo.Var1Low = -3.5
process.demo.Var1High = 3.5


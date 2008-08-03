# The following comments couldn't be translated into the new config version:

#
# this is a sample configuration file for the Tag and Probe Analyser used
# for the CSA08 exercise May-June 2008 --- Nikolaos Rompotis 04.06.08 ---
# 
# ... 15 July, 2008: Kalanand Mishra, Fermilab imported the file into 
#    generic TagAndProbe package for the purpose of comparing the new 
#    tag-and-probe method with the earlier method for electron.
#

#
#keep the logging output to a nice level

import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")
import EgammaAnalysis.ElectronIDProducers.electronIdCutBased_cfi
process.electronIdRobust = EgammaAnalysis.ElectronIDProducers.electronIdCutBased_cfi.eidCutBased.clone()
import EgammaAnalysis.ElectronIDProducers.electronIdCutBased_cfi
process.electronIdLoose = EgammaAnalysis.ElectronIDProducers.electronIdCutBased_cfi.eidCutBased.clone()
import EgammaAnalysis.ElectronIDProducers.electronIdCutBased_cfi
process.electronIdTight = EgammaAnalysis.ElectronIDProducers.electronIdCutBased_cfi.eidCutBased.clone()
import EgammaAnalysis.ElectronIDProducers.electronIdCutBased_cfi
process.electronIdTightRobust = EgammaAnalysis.ElectronIDProducers.electronIdCutBased_cfi.eidCutBased.clone()
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring()
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1)
)
process.MessageLogger = cms.Service("MessageLogger")

process.analysis = cms.EDFilter("TagProbeAnalyserCSA08",
    SCCollectionHybrid = cms.InputTag("correctedHybridSuperClusters"),
    electronIDAssocProducerTightRobust = cms.InputTag("electronIdTightRobust"),
    HLTFilterType = cms.string('hltL1NonIsoHLTNonIsoSingleElectronEt15TrackIsolFilter'),
    EndcapMaxEta = cms.untracked.double(2.5),
    ElectronCollection = cms.InputTag("pixelMatchGsfElectrons"),
    ElectronLabel = cms.string(' '),
    outputfile = cms.untracked.string('./tp.root'),
    # Track Collections
    CtfTrackCollection = cms.InputTag("generalTracks"),
    IsoConeMaxDR = cms.untracked.double(0.6),
    electronIDAssocProducerRobust = cms.InputTag("electronIdRobust"),
    # TIP cutoff: the default value for Claire's 167 version was 0.06
    TagTIPCut = cms.double(999.0),
    TagProbeMassMin = cms.untracked.double(85.0),
    TrackInIsoConeMinPt = cms.untracked.double(1.5),
    RecoEleSeedBCMaxDE = cms.untracked.double(0.0001),
    #  include "EgammaAnalysis/EgammaEfficiencyAlgos/test/CutValues_85_95_LooseIso.cfi" 
    #  these were originally in the file above
    BarrelMaxEta = cms.untracked.double(1.4442),
    IsoConeMinDR = cms.untracked.double(0.02),
    #
    HLTCollection = cms.InputTag("hltTriggerSummaryAOD"),
    GsfTrackMinInnerPt = cms.untracked.double(10.0),
    TagProbeMassMax = cms.untracked.double(95.0),
    electronIDAssocProducerTight = cms.InputTag("electronIdTight"),
    #
    #
    useTriggerInfo = cms.untracked.bool(True),
    TagElectronMinEt = cms.untracked.double(15.0),
    IsoMaxSumPt = cms.untracked.double(0.02),
    ProbeRecoEleSCMaxDE = cms.untracked.double(0.001),
    MCCollection = cms.InputTag("source"),
    useTagIsolation = cms.untracked.bool(True),
    SCCollectionIsland = cms.InputTag("correctedEndcapSuperClustersWithPreshower"),
    ProbeHLTObjMaxDR = cms.untracked.double(0.1),
    EndcapMinEta = cms.untracked.double(1.56),
    ProbeSCMinEt = cms.untracked.double(20.0),
    electronIDAssocProducerLoose = cms.InputTag("electronIdLoose")
)

process.eca = cms.EDAnalyzer("EventContentAnalyzer")

process.runEleID = cms.Sequence(process.electronIdRobust+process.electronIdLoose+process.electronIdTight+process.electronIdTightRobust)
process.main = cms.Sequence(process.runEleID*process.analysis)
process.p = cms.Path(process.main)
process.PoolSource.fileNames = ['file:/uscms_data/d1/kalanand/CSA08-Zee/165E4BA2-AC2B-DD11-8163-001A644EB7CE.root', 'file:/uscms_data/d1/kalanand/CSA08-Zee/2E9A31CC-AC2B-DD11-931E-001A644EB264.root', 'file:/uscms_data/d1/kalanand/CSA08-Zee/5C0A6698-A42B-DD11-B7C5-00145EED0908.root', 'file:/uscms_data/d1/kalanand/CSA08-Zee/8A71884B-992B-DD11-8922-001A644EB2CA.root', 'file:/uscms_data/d1/kalanand/CSA08-Zee/90CB8C85-A02B-DD11-A36F-001A6434EF14.root', 
    'file:/uscms_data/d1/kalanand/CSA08-Zee/A29C770D-AA2B-DD11-9F2F-00145EED0788.root']


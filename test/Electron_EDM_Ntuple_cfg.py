# The following comments couldn't be translated into the new config version:

# keep the logging output to a nice level

import FWCore.ParameterSet.Config as cms

process = cms.Process("ElectronEff")
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

process.load("MuonAnalysis.TagAndProbe.tag_probe_electron_cfi")

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring()
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)
process.MessageLogger = cms.Service("MessageLogger",
    destinations = cms.untracked.vstring('cout', 
        'cerr')
)

process.TPEdm = cms.EDFilter("TagProbeEDMNtuple",
    allProbeTruthMatchMapTags = cms.untracked.VInputTag(cms.InputTag("SuperClustersMatch"), cms.InputTag("GsfElectronsMatch"), cms.InputTag("IsolationMatch"), cms.InputTag("IdMatch")),
    checkExactOverlap = cms.untracked.bool(False),
    triggerDelRMatch = cms.untracked.double(0.3),
    hltL1Tag = cms.untracked.InputTag("hltL1IsoSingleElectronTrackIsolFilter"),
    allProbeCandTags = cms.untracked.VInputTag(cms.InputTag("theSuperClusters"), cms.InputTag("theGsfElectrons"), cms.InputTag("theIsolation"), cms.InputTag("theId")),
    # Trigger & Trigger Matching tags
    triggerEventTag = cms.untracked.InputTag("hltTriggerSummaryRAW"),
    passProbeTruthMatchMapTags = cms.untracked.VInputTag(cms.InputTag("GsfElectronsMatch"), cms.InputTag("IsolationMatch"), cms.InputTag("IdMatch"), cms.InputTag("HLTMatch")),
    # Tag & Probe Electron Candidate Collections
    tagCandTags = cms.untracked.VInputTag(cms.InputTag("theHLT"), cms.InputTag("theHLT"), cms.InputTag("theHLT"), cms.InputTag("theHLT")),
    hltTag = cms.untracked.InputTag("hltL1IsoSingleElectronPixelMatchFilter"),
    # Tag & Probe Muon Association Map 
    tagProbeMapTags = cms.untracked.VInputTag(cms.InputTag("tpMapSuperClusters"), cms.InputTag("tpMapGsfElectrons"), cms.InputTag("tpMapIsolation"), cms.InputTag("tpMapId")),
    # Type of tag-probe candidates, use "Muon" or "Electron"
    # For the moment this only affects the kind of particle
    # used for storing MC truth information.
    tagProbeType = cms.untracked.string('Electron'),
    # Truth Map Tags
    tagTruthMatchMapTags = cms.untracked.VInputTag(cms.InputTag("HLTMatch"), cms.InputTag("HLTMatch"), cms.InputTag("HLTMatch"), cms.InputTag("HLTMatch")),
    # Store some generic information about the event
    # in case we want it
    mcParticles = cms.untracked.vint32(23, 11, 22),
    trackTags = cms.untracked.VInputTag(cms.InputTag("generalTracks")),
    # Pass Probe Electron Candidate Collections
    passProbeCandTags = cms.untracked.VInputTag(cms.InputTag("theGsfElectrons"), cms.InputTag("theIsolation"), cms.InputTag("theId"), cms.InputTag("theHLT")),
    verticesTag = cms.untracked.InputTag("offlinePrimaryVertices"),
    mcParents = cms.untracked.vint32(0, 0, 0)
)

process.outpath = cms.OutputModule("PoolOutputModule",
    outputCommands = cms.untracked.vstring('drop *', 
        'keep *_TPEdm_*_*'),
    #     untracked string fileName = "output/CSA08_electron_tp_edm_ntuple_Reco.root"
    fileName = cms.untracked.string('output/CSA08_electron_tp_edm_ntuple_reReco.root')
)

process.p1 = cms.Path(process.lepton_cands+process.TPEdm)
process.the_end = cms.EndPath(process.outpath)
process.PoolSource.fileNames = ['file:/uscms_data/d1/kalanand/CSA08-Zee/165E4BA2-AC2B-DD11-8163-001A644EB7CE.root', 'file:/uscms_data/d1/kalanand/CSA08-Zee/2E9A31CC-AC2B-DD11-931E-001A644EB264.root', 'file:/uscms_data/d1/kalanand/CSA08-Zee/5C0A6698-A42B-DD11-B7C5-00145EED0908.root', 'file:/uscms_data/d1/kalanand/CSA08-Zee/8A71884B-992B-DD11-8922-001A644EB2CA.root', 'file:/uscms_data/d1/kalanand/CSA08-Zee/90CB8C85-A02B-DD11-A36F-001A6434EF14.root', 
    'file:/uscms_data/d1/kalanand/CSA08-Zee/A29C770D-AA2B-DD11-9F2F-00145EED0788.root']


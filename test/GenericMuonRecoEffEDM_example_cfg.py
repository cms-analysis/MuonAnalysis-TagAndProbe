# The following comments couldn't be translated into the new config version:

# keep the logging output to a nice level

import FWCore.ParameterSet.Config as cms

process = cms.Process("MuonEff")
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

process.load("MuonAnalysis.TagAndProbe.tag_probe_generic_cfi")

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
    allProbeTruthMatchMapTags = cms.untracked.VInputTag(cms.InputTag("allProbeMuonMatch")),
    triggerDelPtRelMatch = cms.untracked.double(0.15),
    triggerDelRMatch = cms.untracked.double(0.15),
    hltL1Tag = cms.untracked.InputTag("hltSingleMuIsoLevel1Seed"),
    allProbeCandTags = cms.untracked.VInputTag(cms.InputTag("probeCands")),
    # Trigger & Trigger Matching tags
    triggerEventTag = cms.untracked.InputTag("hltTriggerSummaryRAW"),
    passProbeTruthMatchMapTags = cms.untracked.VInputTag(cms.InputTag("passProbeMuonMatch")),
    # Tag & Probe Muon Candidate Collections
    tagCandTags = cms.untracked.VInputTag(cms.InputTag("tagCands")),
    hltTag = cms.untracked.InputTag("hltSingleMuIsoL3IsoFiltered"),
    # Tag & Probe Muon Association Map ... can have many, for
    # as many kinds of probe/pass-probe combination as desired.
    # (Then the above muon and cand arrays should be filled appropriately)
    tagProbeMapTags = cms.untracked.VInputTag(cms.InputTag("muonTagProbeMap")),
    # Type of tag-probe candidates, use "Muon" or "Electron"
    # For the moment this only affects the kind of particle
    # used for storing MC truth information.
    tagProbeType = cms.untracked.string('Muon'),
    # Truth Map Tags
    tagTruthMatchMapTags = cms.untracked.VInputTag(cms.InputTag("tagMuonMatch")),
    # Store some generic information about the event
    # in case we want it
    mcParticles = cms.untracked.vint32(23, 13, 22),
    passProbeCandTags = cms.untracked.VInputTag(cms.InputTag("tkStaMatched")),
    mcParents = cms.untracked.vint32(0, 0, 0)
)

process.outpath = cms.OutputModule("PoolOutputModule",
    outputCommands = cms.untracked.vstring('drop *', 
        'keep *_TPEdm_*_*'),
    fileName = cms.untracked.string('output/generic_tp_edm_ntuple_cmssw207_1.root')
)

process.p1 = cms.Path(process.muon_cands+process.TPEdm)
process.the_end = cms.EndPath(process.outpath)
process.PoolSource.fileNames = ['/store/mc/CSA08/Zmumu/GEN-SIM-RECO/1PB_V2_RECO_v1/0028/10375028-DD24-DD11-8170-001617E30F4C.root', '/store/mc/CSA08/Zmumu/GEN-SIM-RECO/1PB_V2_RECO_v1/0028/12CA1A34-DB24-DD11-AA5F-001D09F2426D.root', '/store/mc/CSA08/Zmumu/GEN-SIM-RECO/1PB_V2_RECO_v1/0028/F037A61F-DC24-DD11-B3F4-000423D98844.root', '/store/mc/CSA08/Zmumu/GEN-SIM-RECO/1PB_V2_RECO_v1/0029/3001CA6B-E024-DD11-830E-001D09F292D1.root', '/store/mc/CSA08/Zmumu/GEN-SIM-RECO/1PB_V2_RECO_v1/0029/D06BC70D-E124-DD11-BB0A-001D09F25041.root', 
    '/store/mc/CSA08/Zmumu/GEN-SIM-RECO/1PB_V2_RECO_v1/0029/D6B9828D-DE24-DD11-AE72-001D09F2432B.root']


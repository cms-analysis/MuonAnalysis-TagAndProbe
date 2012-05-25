import FWCore.ParameterSet.Config as cms

from MuonAnalysis.TagAndProbe.mvaIsoVariables_cfi import mvaIsoVariables

from CommonTools.ParticleFlow.pfNoPileUp_cff import *
# default pfNoPU reads from somewhere else. don't want to mess up with it, so make copy
pfPileUpMVAIso = pfPileUp.clone(PFCandidates = 'particleFlow')
pfNoPileUpMVAIso = pfNoPileUp.clone(topCollection = 'pfPileUpMVAIso', bottomCollection = 'particleFlow')

## We keep PfNoPU for use with radialIso
mvaIsoVariablesSeq = cms.Sequence(
    mvaIsoVariables +
    pfPileUpMVAIso * pfNoPileUpMVAIso 
)

MVAIsoVariablesPlain = cms.PSet(
    ChargedIso_DR0p0To0p1 = cms.InputTag("mvaIsoVariables", "ChargedIsoDR0p0To0p1"),
    ChargedIso_DR0p1To0p2 = cms.InputTag("mvaIsoVariables", "ChargedIsoDR0p1To0p2"),
    ChargedIso_DR0p2To0p3 = cms.InputTag("mvaIsoVariables", "ChargedIsoDR0p2To0p3"),
    ChargedIso_DR0p3To0p4 = cms.InputTag("mvaIsoVariables", "ChargedIsoDR0p3To0p4"),
    ChargedIso_DR0p4To0p5 = cms.InputTag("mvaIsoVariables", "ChargedIsoDR0p4To0p5"),
    GammaIso_DR0p0To0p1 = cms.InputTag("mvaIsoVariables", "GammaIsoDR0p0To0p1"),
    GammaIso_DR0p1To0p2 = cms.InputTag("mvaIsoVariables", "GammaIsoDR0p1To0p2"),
    GammaIso_DR0p2To0p3 = cms.InputTag("mvaIsoVariables", "GammaIsoDR0p2To0p3"),
    GammaIso_DR0p3To0p4 = cms.InputTag("mvaIsoVariables", "GammaIsoDR0p3To0p4"),
    GammaIso_DR0p4To0p5 = cms.InputTag("mvaIsoVariables", "GammaIsoDR0p4To0p5"),
    NeutralHadronIso_DR0p0To0p1 = cms.InputTag("mvaIsoVariables", "NeutralHadronIsoDR0p0To0p1"),
    NeutralHadronIso_DR0p1To0p2 = cms.InputTag("mvaIsoVariables", "NeutralHadronIsoDR0p1To0p2"),
    NeutralHadronIso_DR0p2To0p3 = cms.InputTag("mvaIsoVariables", "NeutralHadronIsoDR0p2To0p3"),
    NeutralHadronIso_DR0p3To0p4 = cms.InputTag("mvaIsoVariables", "NeutralHadronIsoDR0p3To0p4"),
    NeutralHadronIso_DR0p4To0p5 = cms.InputTag("mvaIsoVariables", "NeutralHadronIsoDR0p4To0p5"),
    PFCharged = cms.InputTag("mvaIsoVariables", "PFCharged"),
    PFNeutral = cms.InputTag("mvaIsoVariables", "PFNeutral"),
    PFPhotons = cms.InputTag("mvaIsoVariables", "PFPhotons"),
    SumDeltaR  = cms.InputTag("mvaIsoVariables", "SumDeltaR"),
    DeltaRMean = cms.InputTag("mvaIsoVariables", "DeltaRMean"),
    Density    = cms.InputTag("mvaIsoVariables", "Density"),
)

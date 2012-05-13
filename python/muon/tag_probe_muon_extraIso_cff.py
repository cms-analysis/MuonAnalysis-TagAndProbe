import FWCore.ParameterSet.Config as cms

ExtraIsolationVariables = cms.PSet(
    kt6RhoAll         = cms.InputTag("computeCorrectedIso", "Rho"),
    ecalIsoRhoCorr    = cms.InputTag("computeCorrectedIso", "ecalIsoRhoCorrected"),
    hcalIsoRhoCorr    = cms.InputTag("computeCorrectedIso", "hcalIsoRhoCorrected"),
    combRelIsoRhoCorr = cms.InputTag("computeCorrectedIso", "combRelIsoRhoCorrected"),
)

from CommonTools.ParticleFlow.pfNoPileUp_cff import *
# default pfNoPU reads from somewhere else. don't want to mess up with it, so make copy
pfPileUpPFIso = pfPileUp.clone(PFCandidates = 'particleFlow')
pfNoPileUpPFIso = pfNoPileUp.clone(topCollection = 'pfPileUpPFIso', bottomCollection = 'particleFlow')


muonPFIsoChHad04 = cms.EDProducer("MuonPFIsoSingleTypeMapProd",
    muonLabel   = cms.InputTag("probeMuons"),
    pfLabel     = cms.InputTag("pfNoPileUpPFIso"),
    pfSelection = cms.string("charge != 0 && abs(pdgId) == 211"), # charged hadrons
    deltaR      = cms.double(0.4),
)
muonPFIsoNHad04pt0  = muonPFIsoChHad04.clone(pfSelection = "charge == 0 && abs(pdgId) == 130") # neutral hadrons
muonPFIsoNHad04pt05 = muonPFIsoChHad04.clone(pfSelection = "charge == 0 && abs(pdgId) == 130 && pt > 0.5") # 
muonPFIsoNHad04pt1  = muonPFIsoChHad04.clone(pfSelection = "charge == 0 && abs(pdgId) == 130 && pt > 1")   # 
muonPFIsoPhot04pt0  = muonPFIsoChHad04.clone(pfSelection = "charge == 0 && abs(pdgId) == 22") # neutral hadrons
muonPFIsoPhot04pt05 = muonPFIsoChHad04.clone(pfSelection = "charge == 0 && abs(pdgId) == 22 && pt > 0.5") # 
muonPFIsoPhot04pt1  = muonPFIsoChHad04.clone(pfSelection = "charge == 0 && abs(pdgId) == 22 && pt > 1")   # 

muonPFIsoChHad03    = muonPFIsoChHad04.clone(deltaR = 0.3)
muonPFIsoNHad03pt0  = muonPFIsoNHad04pt0.clone(deltaR = 0.3)
muonPFIsoNHad03pt05 = muonPFIsoNHad04pt05.clone(deltaR = 0.3)
muonPFIsoNHad03pt1  = muonPFIsoNHad04pt1.clone(deltaR = 0.3)
muonPFIsoPhot03pt0  = muonPFIsoPhot04pt0.clone(deltaR = 0.3)
muonPFIsoPhot03pt05 = muonPFIsoPhot04pt05.clone(deltaR = 0.3)
muonPFIsoPhot03pt1  = muonPFIsoPhot04pt1.clone(deltaR = 0.3)

muonPFIsoChHad04PU = muonPFIsoChHad04.clone(pfLabel = "pfPileUpPFIso")
muonPFIsoChHad03PU = muonPFIsoChHad03.clone(pfLabel = "pfPileUpPFIso")

muonPFIsoSequence = cms.Sequence(
    pfPileUpPFIso * pfNoPileUpPFIso * (
        muonPFIsoChHad03    +
        muonPFIsoChHad04    +
        muonPFIsoChHad03PU  +
        muonPFIsoChHad04PU  +
        muonPFIsoNHad03pt0  +
        muonPFIsoNHad03pt05 +
        muonPFIsoNHad03pt1  +
        muonPFIsoNHad04pt0  +
        muonPFIsoNHad04pt05 +
        muonPFIsoNHad04pt1  +
        muonPFIsoPhot03pt0  +
        muonPFIsoPhot03pt05 +
        muonPFIsoPhot03pt1  +
        muonPFIsoPhot04pt0  +
        muonPFIsoPhot04pt05 +
        muonPFIsoPhot04pt1  
    )
)

MuonPFIsoVariables = cms.PSet(
   chargedHadIso03  = cms.InputTag("muonPFIsoChHad03"),    
   chargedHadIso04  = cms.InputTag("muonPFIsoChHad04"),    
   puIso03  = cms.InputTag("muonPFIsoChHad03PU"),  
   puIso04  = cms.InputTag("muonPFIsoChHad04PU"), 
   neutralHadIso03Loose  = cms.InputTag("muonPFIsoNHad03pt0"),  
   neutralHadIso03       = cms.InputTag("muonPFIsoNHad03pt05"), 
   neutralHadIso03Tight  = cms.InputTag("muonPFIsoNHad03pt1"),  
   neutralHadIso04Loose  = cms.InputTag("muonPFIsoNHad04pt0"),  
   neutralHadIso04       = cms.InputTag("muonPFIsoNHad04pt05"), 
   neutralHadIso04Tight  = cms.InputTag("muonPFIsoNHad04pt1"),  
   photonIso03Loose  = cms.InputTag("muonPFIsoPhot03pt0"),  
   photonIso03       = cms.InputTag("muonPFIsoPhot03pt05"), 
   photonIso03Tight  = cms.InputTag("muonPFIsoPhot03pt1"),  
   photonIso04Loose  = cms.InputTag("muonPFIsoPhot04pt0"),  
   photonIso04       = cms.InputTag("muonPFIsoPhot04pt05"), 
   photonIso04Tight  = cms.InputTag("muonPFIsoPhot04pt1"),  
)



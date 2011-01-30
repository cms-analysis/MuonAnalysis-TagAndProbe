import FWCore.ParameterSet.Config as cms

from MuonAnalysis.TagAndProbe.nearbyMuonsInfo_cfi import nearbyMuonsInfo as tagProbeSeparation

#########################################################################################
##        Object counting modules                                                      ##
#########################################################################################

nverticesModule = cms.EDProducer("VertexMultiplicityCounter", 
    probes = cms.InputTag("tagMuons"),
    objects = cms.InputTag("offlinePrimaryVertices"),
    objectSelection = cms.string("!isFake && ndof > 4 && abs(z) <= 25 && position.Rho <= 2"),
)

njets15Module = cms.EDProducer("CandCleanedMultiplicityCounter", 
    pairs   = cms.InputTag("tpPairs"),
    objects = cms.InputTag("ak5PFJets"),
    objectSelection = cms.string("abs(eta) < 5 && pt > 15"), 
    minTagObjDR   = cms.double(0.3),
    minProbeObjDR = cms.double(0.3),
)

njets30Module = cms.EDProducer("CandCleanedMultiplicityCounter", 
    pairs   = cms.InputTag("tpPairs"),
    objects = cms.InputTag("ak5PFJets"),
    objectSelection = cms.string("abs(eta) < 5 && pt > 30"), 
    minTagObjDR   = cms.double(0.3),
    minProbeObjDR = cms.double(0.3),
)

#########################################################################################
##        Tracking-related modules                                                     ##
#########################################################################################

## Now I have to define the passing probes for tracking
## first remove low pt tracks which will not make muons anyway
pCutTracks = cms.EDFilter("TrackSelector", 
    src = cms.InputTag("generalTracks"),      
    cut = cms.string("pt > 2 || (abs(eta) > 1 && p > 2)"),
)
tkTracks = cms.EDProducer("ConcreteChargedCandidateProducer", 
    src = cms.InputTag("pCutTracks"),
    particleType = cms.string("mu+"),
)

## Filter out the J/Psi's, to compute fake matching rate
tkTracksNoJPsi = cms.EDProducer("CandidateResonanceInefficiencyCreator",
    src = cms.InputTag("tkTracks"),
    tags = cms.InputTag("tagMuons"),
    mass    = cms.double(3.096),
    massMin = cms.double(2.85), ## Should cut away
    massMax = cms.double(3.25), ## 99.5% of signal
    onlyBestMatch = cms.bool(False),
    outputMode = cms.string("RefToBaseVector"),
)
tkTracksNoBestJPsi = tkTracksNoJPsi.clone(onlyBestMatch = True)

## Filter out the J/Psi's, to compute fake matching rate
tkTracksNoZ = cms.EDProducer("CandidateResonanceInefficiencyCreator",
    src = cms.InputTag("tkTracks"),
    tags = cms.InputTag("tagMuons"),
    mass    = cms.double(91.2),
    massMin = cms.double(40),  ## Should cut away most
    massMax = cms.double(200), ## of the signal
    onlyBestMatch = cms.bool(True),
    outputMode = cms.string("RefToBaseVector"),
)

staToTkMatch = cms.EDProducer("MatcherUsingTracks",
    src     = cms.InputTag("probeMuonsSta"),
    matched = cms.InputTag("tkTracks"),  
    algorithm = cms.string("byDirectComparison"), 
    srcTrack     = cms.string("muon"),    srcState = cms.string("atVertex"), 
    matchedTrack = cms.string("tracker"), matchedState = cms.string("atVertex"),
    maxDeltaR        = cms.double(1.),   # large range in DR (we can tighten it later)
    maxDeltaEta      = cms.double(0.4),  # small in eta, which is more precise
    maxDeltaLocalPos = cms.double(100),
    maxDeltaPtRel    = cms.double(5),   # |pt(sta) - pt(tk)|/pt(tk)
    sortBy           = cms.string("deltaR"),
    requireSameCharge = cms.bool(True),
)
staToTkMatchNoJPsi = staToTkMatch.clone(matched = 'tkTracksNoJPsi')
staToTkMatchNoBestJPsi = staToTkMatch.clone(matched = 'tkTracksNoBestJPsi')
staToTkMatchNoZ = staToTkMatch.clone(matched = 'tkTracksNoZ')

staToTkMatchSequenceJPsi = cms.Sequence(
    pCutTracks + 
    tkTracks           * staToTkMatch           +
    tkTracksNoJPsi     * staToTkMatchNoJPsi     +
    tkTracksNoBestJPsi * staToTkMatchNoBestJPsi 
)
staToTkMatchSequenceZ = cms.Sequence(
    pCutTracks +
    tkTracks    * staToTkMatch    +
    tkTracksNoZ * staToTkMatchNoZ     
)

#########################################################################################
##        Other modules                                                                ##
#########################################################################################

muonDxyPVdzmin = cms.EDProducer("MuonDxyPVdzmin",
    probes = cms.InputTag("probeMuons"),
)

muonsPassingPF = cms.EDProducer("MuonsPassingPF",
    muons = cms.InputTag("probeMuons"),
    pf    = cms.InputTag("particleFlow"),
    matchByReference = cms.bool(False), # set to true only if your probeMuons are a subset by reference of the "muons" collection
)

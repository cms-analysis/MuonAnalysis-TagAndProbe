import FWCore.ParameterSet.Config as cms

KinematicVariables = cms.PSet(
    pt  = cms.string("pt"),
    p   = cms.string("p"),
    eta = cms.string("eta"),
    phi = cms.string("phi"),
    abseta = cms.string("abs(eta)"),
    charge = cms.string("charge")
)
IsolationVariables = cms.PSet(
    tkIso  = cms.string("isolationR03.sumPt"),
    ecalIso = cms.string("isolationR03.emEt"),
    hcalIso = cms.string("isolationR03.hadEt"),
    combRelIso = cms.string("(isolationR03.emEt + isolationR03.hadEt + isolationR03.sumPt)/pt"),
    neutralHadIso = cms.string("neutralHadronIso"),
    chargedHadIso = cms.string("chargedHadronIso"),
    photonIso = cms.string("photonIso"),
    combRelIsoP = cms.string("(neutralHadronIso + chargedHadronIso + photonIso)/pt"),
)
MuonIDVariables = cms.PSet(
    caloCompatibility = cms.string("? isCaloCompatibilityValid ? caloCompatibility : -1"),
    numberOfMatches   = cms.string("? isMatchesValid ? numberOfMatches : -1"),
    numberOfMatchedStations = cms.string("? isMatchesValid ? numberOfMatchedStations : -1"),
)
TrackQualityVariables = cms.PSet(
    dB          = cms.string("dB"),
    tkValidHits = cms.string("? track.isNull ? 0 : track.numberOfValidHits"),
    tkValidPixelHits = cms.string("? track.isNull ? 0 : track.hitPattern.numberOfValidPixelHits"),
    tkPixelLay  = cms.string("? track.isNull ? 0 : track.hitPattern.pixelLayersWithMeasurement"),
    tkExpHitIn  = cms.string("? track.isNull ? 0 : track.trackerExpectedHitsInner.numberOfLostHits"),
    tkExpHitOut = cms.string("? track.isNull ? 0 : track.trackerExpectedHitsOuter.numberOfLostHits"),
    tkHitFract  = cms.string("? track.isNull ? 0 : track.hitPattern.numberOfValidHits/(track.hitPattern.numberOfValidHits+track.hitPattern.numberOfLostHits+track.trackerExpectedHitsInner.numberOfLostHits+track.trackerExpectedHitsOuter.numberOfLostHits)"),
    tkChi2 = cms.string("? track.isNull ? -1 : track.normalizedChi2"),
    tkPtError = cms.string("? track.isNull ? -1 : track.ptError"),
    tkSigmaPtOverPt = cms.string("? track.isNull ? -1 : track.ptError/track.pt"),
)
GlobalTrackQualityVariables = cms.PSet(
    glbChi2 = cms.string("? globalTrack.isNull ? -1 : globalTrack.normalizedChi2"),
    glbPtError = cms.string("? globalTrack.isNull ? -1 : globalTrack.ptError"),
    glbSigmaPtOverPt = cms.string("? globalTrack.isNull ? -1 : globalTrack.ptError/globalTrack.pt"),
)
L1Variables = cms.PSet(
    l1pt = cms.string("? userCand('muonL1Info').isNull ? 0 : userCand('muonL1Info').pt"),
    l1q  = cms.string("userInt('muonL1Info:quality')"),
    l1dr = cms.string("userFloat('muonL1Info:deltaR')"),
)
L2Variables = cms.PSet(
    l2pt  = cms.string("? triggerObjectMatchesByCollection('hltL2MuonCandidates').empty() ? 0 : triggerObjectMatchesByCollection('hltL2MuonCandidates').at(0).pt"),
    l2eta = cms.string("? triggerObjectMatchesByCollection('hltL2MuonCandidates').empty() ? 0 : triggerObjectMatchesByCollection('hltL2MuonCandidates').at(0).eta"),
    l2dr  = cms.string("? triggerObjectMatchesByCollection('hltL2MuonCandidates').empty() ? 999 : "+
                      " deltaR( eta, phi, " +
                      "         triggerObjectMatchesByCollection('hltL2MuonCandidates').at(0).eta, "+
                      "         triggerObjectMatchesByCollection('hltL2MuonCandidates').at(0).phi ) ")
)
L3Variables = cms.PSet(
    l3pt = cms.string("? triggerObjectMatchesByCollection('hltL3MuonCandidates').empty() ? 0 : triggerObjectMatchesByCollection('hltL3MuonCandidates').at(0).pt"),
    l3dr = cms.string("? triggerObjectMatchesByCollection('hltL3MuonCandidates').empty() ? 999 : "+
                      " deltaR( eta, phi, " +
                      "         triggerObjectMatchesByCollection('hltL3MuonCandidates').at(0).eta, "+
                      "         triggerObjectMatchesByCollection('hltL3MuonCandidates').at(0).phi ) ")
)
TriggerVariables = cms.PSet(L1Variables, L2Variables, L3Variables)
AllVariables = cms.PSet(KinematicVariables, IsolationVariables, MuonIDVariables, TrackQualityVariables, GlobalTrackQualityVariables, L1Variables, L2Variables, L3Variables)

TrackQualityFlags = cms.PSet(
    Track_HP  = cms.string("? track.isNonnull ? track.quality('highPurity') : 0"),
)
MuonIDFlags = cms.PSet(
    Calo   = cms.string("isCaloMuon"),
    Glb    = cms.string("isGlobalMuon"),
    GlbPT  = cms.string("muonID('GlobalMuonPromptTight')"),
    TM     = cms.string("isTrackerMuon"),
    TMA    = cms.string("muonID('TrackerMuonArbitrated')"),
    TMLSAT = cms.string("muonID('TMLastStationAngTight')"),
    TMOSL  = cms.string("muonID('TMOneStationLoose')"),
    TMOST  = cms.string("muonID('TMOneStationTight')"),
    VBTF   = cms.string("numberOfMatchedStations > 1 && muonID('GlobalMuonPromptTight') && abs(dB) < 0.2 && "+
                        "track.numberOfValidHits > 10 && track.hitPattern.numberOfValidPixelHits > 0"),
)



HighPtTriggerFlags = cms.PSet(

   # legacy
   Mu9       = cms.string("!triggerObjectMatchesByPath('HLT_Mu9').empty()"),
   DoubleMu3 = cms.string("!triggerObjectMatchesByPath('HLT_DoubleMu3_v*').empty()"),
   # current or new up to 3E33 menu
   Mu15      = cms.string("!triggerObjectMatchesByPath('HLT_Mu15_v*',1,0).empty()"),
   Mu20      = cms.string("!triggerObjectMatchesByPath('HLT_Mu20_v*',1,0).empty()"),
   Mu24      = cms.string("!triggerObjectMatchesByPath('HLT_Mu24_v*',1,0).empty()"),
   Mu30      = cms.string("!triggerObjectMatchesByPath('HLT_Mu30_v*',1,0).empty()"),
   Mu40      = cms.string("!triggerObjectMatchesByPath('HLT_Mu40_v*',1,0).empty()"),
   IsoMu15   = cms.string("!triggerObjectMatchesByPath('HLT_IsoMu15_v*',1,0).empty()"),
   IsoMu17   = cms.string("!triggerObjectMatchesByPath('HLT_IsoMu17_v*',1,0).empty()"),
   IsoMu24_eta2p1   = cms.string("!triggerObjectMatchesByPath('HLT_IsoMu24_eta2p1_v*',1,0).empty()"),
   IsoMu30_eta2p1   = cms.string("!triggerObjectMatchesByPath('HLT_IsoMu30_eta2p1_v*',1,0).empty()"),
   IsoMu24   = cms.string("!triggerObjectMatchesByPath('HLT_IsoMu24_v*',1,0).empty()"),
   IsoMu30   = cms.string("!triggerObjectMatchesByPath('HLT_IsoMu30_v*',1,0).empty()"),
   DoubleMu5 = cms.string("!triggerObjectMatchesByPath('HLT_DoubleMu5_v*',1,0).empty()"),
   DoubleMu7 = cms.string("!triggerObjectMatchesByPath('HLT_DoubleMu7_v*',1,0).empty()"),
   DoubleMu13Mu8_Mu13 = cms.string("!triggerObjectMatchesByPath('HLT_Mu13_Mu8_v*',1,0).empty()"),
   DoubleMu13Mu8_Mu8 = cms.string("!triggerObjectMatchesByFilter('hltDiMuonL3PreFiltered8').empty() || !triggerObjectMatchesByFilter('hltDiMuonL3p5PreFiltered8').empty()"),
   DoubleMu17Mu8_Mu17 = cms.string("!triggerObjectMatchesByPath('HLT_Mu17_Mu8_v*',1,0).empty()"),
   DoubleMu17Mu8_Mu8 = cms.string("!triggerObjectMatchesByFilter('hltDiMuonL3PreFiltered8').empty() || !triggerObjectMatchesByFilter('hltDiMuonL3p5PreFiltered8').empty()"),
   DoubleMu17TkMu8_TkMu8 = cms.string("!triggerObjectMatchesByPath('HLT_Mu17_TkMu8_v*',1,0).empty()"),
   DoubleMu17TkMu8_Mu17 = cms.string("!triggerObjectMatchesByFilter('hltL3Mu17FromDiMuonFiltered17').empty()")
   )


LowPtTriggerFlagsPhysics = cms.PSet(
   DoubleMu0_Quarkonium  = cms.string("!triggerObjectMatchesByPath('HLT_DoubleMu0',1,0).empty() || "+
                                      "!triggerObjectMatchesByPath('HLT_DoubleMu0_Quarkonium',1,0).empty() || "+
                                      "!triggerObjectMatchesByPath('HLT_DoubleMu0_Quarkonium_v1',1,0).empty()"),
   DoubleMu3_Jpsi_A      = cms.string("!triggerObjectMatchesByFilter('hltDiMuonL3PreFiltered3Jpsi').empty()"),
   DoubleMu3_Jpsi_B      = cms.string("!triggerObjectMatchesByFilter('hltDoubleMu3JpsiL3Filtered').empty()"),
   DoubleMu3_Jpsi        = cms.string("!triggerObjectMatchesByPath('HLT_DoubleMu3_Jpsi_v*',1,0).empty()"),

)

LowPtTriggerFlagsEfficienciesTag = cms.PSet(
   ## Mu + Track will be added automatically (see below)
   ## Mu + L2Mu
   Mu5_L2Mu0_MU = cms.string("!triggerObjectMatchesByCollection('hltL3MuonCandidates').empty() && "+
                             " triggerObjectMatchesByCollection('hltL3MuonCandidates').at(0).hasFilterLabel('hltMu5L2Mu0L3Filtered5')"),
   Mu5_L2Mu2_Jpsi_MU = cms.string("!triggerObjectMatchesByCollection('hltL3MuonCandidates').empty() && "+
                                  " triggerObjectMatchesByCollection('hltL3MuonCandidates').at(0).hasFilterLabel('hltMu5L2Mu2JpsiTrackMassFiltered')"),
)

LowPtTriggerFlagsEfficienciesProbe = cms.PSet(
   ## Mu + Track will be added automatically (see below)
   ## Mu + L2Mu
   Mu5_L2Mu0_L2 = cms.string("!triggerObjectMatchesByCollection('hltL2MuonCandidates').empty() && "+
                             " triggerObjectMatchesByCollection('hltL2MuonCandidates').at(0).hasFilterLabel('hltDiMuonL2PreFiltered0')"),
   Mu5_L2Mu2_Jpsi_L2 = cms.string("!triggerObjectMatchesByCollection('hltL2MuonCandidates').empty() && "+
                             " triggerObjectMatchesByCollection('hltL2MuonCandidates').at(0).hasFilterLabel('hltMu5L2Mu2JpsiTrackMassFiltered')"),
)
for (ptMu,ptTk) in [ (5,0), (5,2), (3,3), (7,5), (7,7) ]:
   filter = "Mu%dTrack%d" % (ptMu,ptTk) if ptTk > 0 else "Mu%dTrack" % ptMu
   setattr(LowPtTriggerFlagsEfficienciesTag, "Mu%d_Track%d_Jpsi_MU" % (ptMu,ptTk), 
            cms.string("!triggerObjectMatchesByCollection('hltL3MuonCandidates').empty() && "+
                       " triggerObjectMatchesByCollection('hltL3MuonCandidates').at(0).hasFilterLabel('hlt"+filter+"JpsiTrackMassFiltered')"))
   setattr(LowPtTriggerFlagsEfficienciesProbe, "Mu%d_Track%d_Jpsi_TK" % (ptMu,ptTk),
            cms.string("!triggerObjectMatchesByCollection('hltMuTrackJpsiCtfTrackCands').empty() && "+
                       " triggerObjectMatchesByCollection('hltMuTrackJpsiCtfTrackCands').at(0).hasFilterLabel('hlt"+filter+"JpsiTrackMassFiltered')"))

LowPtTriggerFlagsEfficiencies = cms.PSet(LowPtTriggerFlagsEfficienciesTag,LowPtTriggerFlagsEfficienciesProbe)

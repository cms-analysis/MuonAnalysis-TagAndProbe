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
    
    chargedHadIso03 = cms.string("pfIsolationR03().sumChargedHadronPt"),
    puIso03 = cms.string("pfIsolationR03().sumPUPt"),
    neutralHadIso03 = cms.string("pfIsolationR03().sumNeutralHadronEt"),
#    neutralHadIso03HT = cms.string("pfIsolationR03().sumNeutralHadronEtHighThreshold"),
    photonIso03 = cms.string("pfIsolationR03().sumPhotonEt"),
#    photonIso03HT = cms.string("pfIsolationR03().sumPhotonEtHighThreshold"),
    chargedParticleIso03 = cms.string("pfIsolationR03().sumChargedParticlePt"),
    combRelIsoPF03 = cms.string("(pfIsolationR03().sumChargedHadronPt + pfIsolationR03().sumNeutralHadronEt + pfIsolationR03().sumPhotonEt)/pt"),
    combRelIsoPF03dBeta = cms.string("(pfIsolationR03().sumChargedHadronPt + max(pfIsolationR03().sumNeutralHadronEt + pfIsolationR03().sumPhotonEt - pfIsolationR03().sumPUPt/2,0.0))/pt"),
#    combRelIsoPF03HT = cms.string("(pfIsolationR03().sumChargedHadronPt + pfIsolationR03().sumNeutralHadronEtHighThreshold + pfIsolationR03().sumPhotonEtHighThreshold)/pt"),

    chargedHadIso04 = cms.string("pfIsolationR04().sumChargedHadronPt"),
    puIso04 = cms.string("pfIsolationR04().sumPUPt"),
    neutralHadIso04 = cms.string("pfIsolationR04().sumNeutralHadronEt"),
#    neutralHadIso04HT = cms.string("pfIsolationR04().sumNeutralHadronEtHighThreshold"),
    photonIso04 = cms.string("pfIsolationR04().sumPhotonEt"),
#    photonIso04HT = cms.string("pfIsolationR04().sumPhotonEtHighThreshold"),
    chargedParticleIso04 = cms.string("pfIsolationR04().sumChargedParticlePt"),
    combRelIsoPF04 = cms.string("(pfIsolationR04().sumChargedHadronPt + pfIsolationR04().sumNeutralHadronEt + pfIsolationR04().sumPhotonEt)/pt"),
    combRelIsoPF04dBeta = cms.string("(pfIsolationR04().sumChargedHadronPt + max(pfIsolationR04().sumNeutralHadronEt + pfIsolationR04().sumPhotonEt - pfIsolationR04().sumPUPt/2,0.0))/pt"),
#    combRelIsoPF04HT = cms.string("(pfIsolationR04().sumChargedHadronPt + pfIsolationR04().sumNeutralHadronEtHighThreshold + pfIsolationR04().sumPhotonEtHighThreshold)/pt"),

#     neutralHadIso = cms.string("neutralHadronIso"),
#     chargedHadIso = cms.string("chargedHadronIso"),
#     photonIso = cms.string("photonIso"),
#     combRelIsoP = cms.string("(neutralHadronIso + chargedHadronIso + photonIso)/pt"),
)

MuonIDVariables = cms.PSet(
    caloCompatibility = cms.string("? isCaloCompatibilityValid ? caloCompatibility : -1"),
    numberOfMatches   = cms.string("? isMatchesValid ? numberOfMatches : -1"),
    numberOfMatchedStations = cms.string("? isMatchesValid ? numberOfMatchedStations : -1"),
    segmentCompatibility = cms.string("segmentCompatibility"),
)
MuonCaloVariables = cms.PSet(
    hadEnergy   = cms.string("calEnergy.had"),
    emEnergy    = cms.string("calEnergy.em"),
    hadS9Energy = cms.string("calEnergy.hadS9"),
    emS9Energy  = cms.string("calEnergy.emS9"),
)
TrackQualityVariables = cms.PSet(
    # 2D variables
    dB          = cms.string("dB"),
    edB         = cms.string("edB"),
    # 3D variables
    IP = cms.string('abs(dB("PV3D"))'),
    IPError = cms.string('edB("PV3D")'),
    SIP = cms.string('abs(dB("PV3D")/edB("PV3D"))'),
    # Hits and such
    tkValidHits = cms.string("? track.isNull ? 0 : track.numberOfValidHits"),
    tkTrackerLay = cms.string("? track.isNull ? 0 : track.hitPattern.trackerLayersWithMeasurement"),
    tkValidPixelHits = cms.string("? track.isNull ? 0 : track.hitPattern.numberOfValidPixelHits"),
    tkPixelLay  = cms.string("? track.isNull ? 0 : track.hitPattern.pixelLayersWithMeasurement"),
    tkExpHitIn  = cms.string("? track.isNull ? 0 : track.trackerExpectedHitsInner.numberOfLostHits"),
    tkExpHitOut = cms.string("? track.isNull ? 0 : track.trackerExpectedHitsOuter.numberOfLostHits"),
    tkHitFract  = cms.string("? track.isNull ? 0 : track.hitPattern.numberOfValidHits/(track.hitPattern.numberOfValidHits+track.hitPattern.numberOfLostHits+track.trackerExpectedHitsInner.numberOfLostHits+track.trackerExpectedHitsOuter.numberOfLostHits)"),
    tkChi2 = cms.string("? track.isNull ? -1 : track.normalizedChi2"),
    tkPtError = cms.string("? track.isNull ? -1 : track.ptError"),
    tkSigmaPtOverPt = cms.string("? track.isNull ? -1 : track.ptError/track.pt"),
    tkKink = cms.string("combinedQuality.trkKink"),
)
GlobalTrackQualityVariables = cms.PSet(
    glbChi2 = cms.string("? globalTrack.isNull ? -1 : globalTrack.normalizedChi2"),
    glbValidMuHits = cms.string("? globalTrack.isNull ? 0 : globalTrack.hitPattern.numberOfValidMuonHits"),
    glbPtError = cms.string("? globalTrack.isNull ? -1 : globalTrack.ptError"),
    glbSigmaPtOverPt = cms.string("? globalTrack.isNull ? -1 : globalTrack.ptError/globalTrack.pt"),
)
StaOnlyVariables = cms.PSet(
    staQoverP      = cms.string("? outerTrack.isNull() ? 0 : outerTrack.qoverp"),
    staQoverPerror = cms.string("? outerTrack.isNull() ? 0 : outerTrack.qoverpError"),
    staValidStations = cms.string("? outerTrack.isNull() ? -1 : outerTrack.hitPattern.muonStationsWithValidHits()"),
)
L1Variables = cms.PSet(
    l1pt = cms.string("? userCand('muonL1Info').isNull ? 0 : userCand('muonL1Info').pt"),
    l1q  = cms.string("userInt('muonL1Info:quality')"),
    l1dr = cms.string("userFloat('muonL1Info:deltaR')"),
    l1ptByQ = cms.string("? userCand('muonL1Info:ByQ').isNull ? 0 : userCand('muonL1Info:ByQ').pt"),
    l1qByQ  = cms.string("userInt('muonL1Info:qualityByQ')"),
    l1drByQ = cms.string("userFloat('muonL1Info:deltaRByQ')"),
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
AllVariables = cms.PSet(KinematicVariables, IsolationVariables, MuonIDVariables, MuonCaloVariables, TrackQualityVariables, GlobalTrackQualityVariables, L1Variables, L2Variables, L3Variables)

TrackQualityFlags = cms.PSet(
    Track_HP  = cms.string("? track.isNonnull ? track.quality('highPurity') : 0"),
)
MuonIDFlags = cms.PSet(
    Calo   = cms.string("isCaloMuon"),
    Glb    = cms.string("isGlobalMuon"),
    GlbPT  = cms.string("muonID('GlobalMuonPromptTight')"),
    TM     = cms.string("isTrackerMuon"),
    TMA    = cms.string("muonID('TrackerMuonArbitrated')"),
    PF     = cms.string("isPFMuon()"),
    TMLSAT = cms.string("muonID('TMLastStationAngTight')"),
    TMLST  = cms.string("muonID('TMLastStationTight')"),
    TMOSL  = cms.string("muonID('TMOneStationLoose')"),
    TMOST  = cms.string("muonID('TMOneStationTight')"),
    TMOSTQual  = cms.string("muonID('TMOneStationTight') && track.numberOfValidHits > 10 && track.normalizedChi2()<1.8 && track.hitPattern.pixelLayersWithMeasurement>1"),
    VBTF   = cms.string("numberOfMatchedStations > 1 && muonID('GlobalMuonPromptTight') && abs(dB) < 0.2 && "+
                        "track.numberOfValidHits > 10 && track.hitPattern.numberOfValidPixelHits > 0"),
    VBTF_nL8    = cms.string("numberOfMatchedStations > 1 && muonID('GlobalMuonPromptTight') && abs(dB) < 0.2 && "+
                        "track.hitPattern.trackerLayersWithMeasurement > 8 && track.hitPattern.numberOfValidPixelHits > 0"),
    VBTF_nL9    = cms.string("numberOfMatchedStations > 1 && muonID('GlobalMuonPromptTight') && abs(dB) < 0.2 && "+
                        "track.hitPattern.trackerLayersWithMeasurement > 9 && track.hitPattern.numberOfValidPixelHits > 0"),
    Tight2012   = cms.string("isPFMuon && numberOfMatchedStations > 1 && muonID('GlobalMuonPromptTight') && abs(dB) < 0.2 && "+
                        "track.hitPattern.trackerLayersWithMeasurement > 5 && track.hitPattern.numberOfValidPixelHits > 0"),
)

HighPtTriggerFlags = cms.PSet(
   # legacy
   #Mu9       = cms.string("!triggerObjectMatchesByPath('HLT_Mu9').empty()"),
   #DoubleMu3 = cms.string("!triggerObjectMatchesByPath('HLT_DoubleMu3_v*').empty()"),
   # current or new up to 3E33 menu
   #Mu15      = cms.string("!triggerObjectMatchesByPath('HLT_Mu15_v*',1,0).empty()"),
   #Mu20      = cms.string("!triggerObjectMatchesByPath('HLT_Mu20_v*',1,0).empty()"),
   #Mu24      = cms.string("!triggerObjectMatchesByPath('HLT_Mu24_v*',1,0).empty()"),
   #Mu24_eta2p1      = cms.string("!triggerObjectMatchesByPath('HLT_Mu24_eta2p1_v*',1,0).empty()"),
   #Mu30      = cms.string("!triggerObjectMatchesByPath('HLT_Mu30_v*',1,0).empty()"),
   #Mu30_eta2p1      = cms.string("!triggerObjectMatchesByPath('HLT_Mu30_eta2p1_v*',1,0).empty()"),
   Mu40      = cms.string("!triggerObjectMatchesByPath('HLT_Mu40_v*',1,0).empty()"),
   Mu40_eta2p1 = cms.string("!triggerObjectMatchesByPath('HLT_Mu40_eta2p1_v*',1,0).empty()"),
   #IsoMu15   = cms.string("!triggerObjectMatchesByPath('HLT_IsoMu15_v*',1,0).empty()"),
   #IsoMu15_eta2p1 = cms.string("!triggerObjectMatchesByPath('HLT_IsoMu15_eta2p1_v*',1,0).empty()"),
   #IsoMu17   = cms.string("!triggerObjectMatchesByPath('HLT_IsoMu17_v*',1,0).empty()"),
   #IsoMu24   = cms.string("!triggerObjectMatchesByPath('HLT_IsoMu24_v*',1,0).empty()"),
   IsoMu24_eta2p1   = cms.string("!triggerObjectMatchesByPath('HLT_IsoMu24_eta2p1_v*',1,0).empty()"),
   IsoMu30   = cms.string("!triggerObjectMatchesByPath('HLT_IsoMu30_v*',1,0).empty()"),
   #IsoMu30_eta2p1   = cms.string("!triggerObjectMatchesByPath('HLT_IsoMu30_eta2p1_v*',1,0).empty()"),
   #Mu15orMu17orMu20orMu24orMu30orMu40   = cms.string("!triggerObjectMatchesByPath('HLT_Mu15_v*',1,0).empty() || !triggerObjectMatchesByPath('HLT_IsoMu17_v*',1,0).empty() || !triggerObjectMatchesByPath('HLT_Mu20_v*',1,0).empty() || !triggerObjectMatchesByPath('HLT_IsoMu20_v*',1,0).empty() || !triggerObjectMatchesByPath('HLT_Mu24_v*',1,0).empty() || !triggerObjectMatchesByPath('HLT_Mu24_eta2p1_v*',1,0).empty() || !triggerObjectMatchesByPath('HLT_IsoMu24_eta2p1_v*',1,0).empty() || !triggerObjectMatchesByPath('HLT_IsoMu24_v*',1,0).empty() || !triggerObjectMatchesByPath('HLT_IsoMu30_eta2p1_v*',1,0).empty() || !triggerObjectMatchesByPath('HLT_IsoMu30_eta2p1_v*',1,0).empty() || !triggerObjectMatchesByPath('HLT_Mu30_v*',1,0).empty() || !triggerObjectMatchesByPath('HLT_Mu30_eta2p1_v*',1,0).empty() || !triggerObjectMatchesByPath('HLT_Mu40_v*',1,0).empty() || !triggerObjectMatchesByPath('HLT_Mu40_eta2p1_v*',1,0).empty()"),
   #DoubleMu5 = cms.string("!triggerObjectMatchesByPath('HLT_DoubleMu5_v*',1,0).empty()"),
   #DoubleMu7 = cms.string("!triggerObjectMatchesByPath('HLT_DoubleMu7_v*',1,0).empty()"),
   DoubleMu13Mu8_Mu13 = cms.string("!triggerObjectMatchesByPath('HLT_Mu13_Mu8_v*',1,0).empty()"),
   DoubleMu13Mu8_Mu8 = cms.string("!triggerObjectMatchesByFilter('hltDiMuonL3PreFiltered8').empty() || !triggerObjectMatchesByFilter('hltDiMuonL3p5PreFiltered8').empty()"),

   ## Heavily prescaled but still useful   
   Mu17 = cms.string("!triggerObjectMatchesByPath('HLT_Mu8_v*',1,0).empty()"),
   Mu8  = cms.string("!triggerObjectMatchesByPath('HLT_Mu17_v*',1,0).empty()"),

   # 2011 version
   # DoubleMu17Mu8_Mu17 = cms.string("!triggerObjectMatchesByPath('HLT_Mu17_Mu8_v*',1,0).empty()"),
   # DoubleMu17Mu8_Mu8 = cms.string("!triggerObjectMatchesByFilter('hltDiMuonL3PreFiltered8').empty() || !triggerObjectMatchesByFilter('hltDiMuonL3p5PreFiltered8').empty()"),
   # DoubleMu17TkMu8_TkMu8 = cms.string("!triggerObjectMatchesByPath('HLT_Mu17_TkMu8_v*',1,0).empty()"),
   # DoubleMu17TkMu8_Mu17 = cms.string("!triggerObjectMatchesByFilter('hltL3Mu17FromDiMuonFiltered17').empty()")


   # To take care of the presence or not of the dz filter, we are requiring three flags for each of the main un-prescaled DoubleMuon Triggers
   # Mu17 leg, Mu8 leg, and (Mu17 leg && fired path)

   DoubleMu17Mu8_Mu17 = cms.string("!triggerObjectMatchesByPath('HLT_Mu17_Mu8_v*',1,0).empty() && (!triggerObjectMatchesByFilter('hltL3fL1DoubleMu10MuOpenL1f0L2f10L3Filtered17').empty() || !triggerObjectMatchesByFilter('hltL3fL1DoubleMu10MuOpenOR3p5L1f0L2f10L3Filtered17').empty())"),
   DoubleMu17Mu8_Mu8 = cms.string("!triggerObjectMatchesByPath('HLT_Mu17_Mu8_v*',1,0).empty() && (!triggerObjectMatchesByFilter('hltL3pfL1DoubleMu10MuOpenL1f0L2pf0L3PreFiltered8').empty() || !triggerObjectMatchesByFilter('hltL3pfL1DoubleMu10MuOpenOR3p5L1f0L2pf0L3PreFiltered8').empty())"),
   DoubleMu17Mu8_Mu17leg = cms.string("!triggerObjectMatchesByFilter('hltL3fL1DoubleMu10MuOpenL1f0L2f10L3Filtered17').empty() || !triggerObjectMatchesByFilter('hltL3fL1DoubleMu10MuOpenOR3p5L1f0L2f10L3Filtered17').empty()"),
   DoubleMu17Mu8_Mu8leg = cms.string("!triggerObjectMatchesByFilter('hltL3pfL1DoubleMu10MuOpenL1f0L2pf0L3PreFiltered8').empty() || !triggerObjectMatchesByFilter('hltL3pfL1DoubleMu10MuOpenOR3p5L1f0L2pf0L3PreFiltered8').empty()"), 


   DoubleMu17TkMu8_Mu17 = cms.string("!triggerObjectMatchesByPath('HLT_Mu17_TkMu8_v*',1,0).empty() && (!triggerObjectMatchesByFilter('hltL3fL1sMu10MuOpenL1f0L2f10L3Filtered17').empty() || !triggerObjectMatchesByFilter('hltL3fL1DoubleMu10MuOpenOR3p5L1f0L2f10L3Filtered17').empty())"),
   DoubleMu17TkMu8_TkMu8 = cms.string("!triggerObjectMatchesByPath('HLT_Mu17_TkMu8_v*',1,0).empty() && !triggerObjectMatchesByFilter('hltDiMuonGlbFiltered17TrkFiltered8').empty()"),
   DoubleMu17TkMu8_Mu17leg = cms.string("!triggerObjectMatchesByFilter('hltL3fL1sMu10MuOpenL1f0L2f10L3Filtered17').empty() || !triggerObjectMatchesByFilter('hltL3fL1DoubleMu10MuOpenOR3p5L1f0L2f10L3Filtered17').empty()"),
   DoubleMu17TkMu8_TkMu8leg = cms.string("!triggerObjectMatchesByFilter('hltDiMuonGlbFiltered17TrkFiltered8').empty()"),
)
HighPtTriggerFlagsDebug = cms.PSet(
   # --- the ones commented out don't save tags ---
   #L1DoubleMu10MuOpenL1Filtered0 = cms.string("!triggerObjectMatchesByFilter('hltL1DoubleMu10MuOpenL1Filtered0').empty()"),
   #L1DoubleMu10MuOpenOR3p5L1Filtered0 = cms.string("!triggerObjectMatchesByFilter('hltL1DoubleMu10MuOpenOR3p5L1Filtered0').empty()"),
   #L1fL1sDoubleMu10MuOpenL1Filtered0 = cms.string("!triggerObjectMatchesByFilter('hltL1fL1sDoubleMu10MuOpenL1Filtered0').empty()"),
   #L1fL1sDoubleMu10MuOpenOR3p5L1Filtered0  = cms.string("!triggerObjectMatchesByFilter('hltL1fL1sDoubleMu10MuOpenOR3p5L1Filtered0').empty()"),
   #L1fL1sMu16Eta2p1L1Filtered0 = cms.string("!triggerObjectMatchesByFilter('hltL1fL1sMu16Eta2p1L1Filtered0').empty()"),
   #L1fL1sMu16L1Filtered0 = cms.string("!triggerObjectMatchesByFilter('hltL1fL1sMu16L1Filtered0').empty()"),
   L1sMu16 = cms.string("!triggerObjectMatchesByFilter('hltL1sMu16').empty()"),
   L1sMu16Eta2p1 = cms.string("!triggerObjectMatchesByFilter('hltL1sMu16Eta2p1').empty()"),
   L1sL1DoubleMu10MuOpen = cms.string("!triggerObjectMatchesByFilter('hltL1sL1DoubleMu10MuOpen').empty()"),
   L1sL1DoubleMu10MuOpenORDoubleMu103p5 = cms.string("!triggerObjectMatchesByFilter('hltL1sL1DoubleMu10MuOpenORDoubleMu103p5').empty()"),
   L2fL1DoubleMu10MuOpenL1f0L2Filtered10 = cms.string("!triggerObjectMatchesByFilter('hltL2fL1DoubleMu10MuOpenL1f0L2Filtered10').empty()"),
   L2fL1DoubleMu10MuOpenOR3p5L1f0L2Filtered10 = cms.string("!triggerObjectMatchesByFilter('hltL2fL1DoubleMu10MuOpenOR3p5L1f0L2Filtered10').empty()"),
   L2fL1sDoubleMu10MuOpenL1f0L2Filtered10 = cms.string("!triggerObjectMatchesByFilter('hltL2fL1sDoubleMu10MuOpenL1f0L2Filtered10').empty()"),
   L2fL1sDoubleMu10MuOpenOR3p5L1f0L2Filtered10 = cms.string("!triggerObjectMatchesByFilter('hltL2fL1sDoubleMu10MuOpenOR3p5L1f0L2Filtered10').empty()"),
   L2fL1sMu16Eta2p1L1f0L2Filtered16Q = cms.string("!triggerObjectMatchesByFilter('hltL2fL1sMu16Eta2p1L1f0L2Filtered16Q').empty()"),
   L2fL1sMu16L1f0L2Filtered16Q = cms.string("!triggerObjectMatchesByFilter('hltL2fL1sMu16L1f0L2Filtered16Q').empty()"),
   L2pfL1DoubleMu10MuOpenL1f0L2PreFiltered0 = cms.string("!triggerObjectMatchesByFilter('hltL2pfL1DoubleMu10MuOpenL1f0L2PreFiltered0').empty()"),
   L2pfL1DoubleMu10MuOpenOR3p5L1f0L2PreFiltered0 = cms.string("!triggerObjectMatchesByFilter('hltL2pfL1DoubleMu10MuOpenOR3p5L1f0L2PreFiltered0').empty()"),
   L3fL1sMu16Eta2p1L1f0L2f16QL3Filtered24Q = cms.string("!triggerObjectMatchesByFilter('hltL3fL1sMu16Eta2p1L1f0L2f16QL3Filtered24Q').empty()"),
   L3fL1sMu16L1f0L2f16QL3Filtered24Q = cms.string("!triggerObjectMatchesByFilter('hltL3fL1sMu16L1f0L2f16QL3Filtered24Q').empty()"),
)


LowPtTriggerFlagsPhysics = cms.PSet(
   #DoubleMu0_Quarkonium  = cms.string("!triggerObjectMatchesByPath('HLT_DoubleMu0',1,0).empty() || "+
   #                                   "!triggerObjectMatchesByPath('HLT_DoubleMu0_Quarkonium',1,0).empty() || "+
   #                                   "!triggerObjectMatchesByPath('HLT_DoubleMu0_Quarkonium_v1',1,0).empty()"),
   #DoubleMu3_Jpsi_A      = cms.string("!triggerObjectMatchesByFilter('hltDiMuonL3PreFiltered3Jpsi').empty()"),
   #DoubleMu3_Jpsi_B      = cms.string("!triggerObjectMatchesByFilter('hltDoubleMu3JpsiL3Filtered').empty()"),
   #DoubleMu3_Jpsi        = cms.string("!triggerObjectMatchesByPath('HLT_DoubleMu3_Jpsi_v*',1,0).empty()"),
)

LowPtTriggerFlagsEfficienciesTag = cms.PSet(
   ## Mu + Track
   Mu5_Track2_Jpsi_MU = cms.string("!triggerObjectMatchesByCollection('hltL3MuonCandidates').empty() && "+
                                   " triggerObjectMatchesByCollection('hltL3MuonCandidates').at(0).hasFilterLabel('hltMu5Track2JpsiTrackMassFiltered')"),
   Mu5_Track3p5_Jpsi_MU = cms.string("!triggerObjectMatchesByCollection('hltL3MuonCandidates').empty() && "+
                                   " triggerObjectMatchesByCollection('hltL3MuonCandidates').at(0).hasFilterLabel('hltMu5Track3p5JpsiTrackMassFiltered')"),
   Mu7_Track7_Jpsi_MU = cms.string("!triggerObjectMatchesByCollection('hltL3MuonCandidates').empty() && "+
                                   " triggerObjectMatchesByCollection('hltL3MuonCandidates').at(0).hasFilterLabel('hltMu7Track7JpsiTrackMassFiltered')"),
   ## Mu + L2Mu
   Mu5_L2Mu0_MU = cms.string("!triggerObjectMatchesByCollection('hltL3MuonCandidates').empty() && "+
                             " triggerObjectMatchesByCollection('hltL3MuonCandidates').at(0).hasFilterLabel('hltMu5L2Mu0L3Filtered5')"),
   Mu5_L2Mu2_Jpsi_MU = cms.string("!triggerObjectMatchesByCollection('hltL3MuonCandidates').empty() && "+
                                  " triggerObjectMatchesByCollection('hltL3MuonCandidates').at(0).hasFilterLabel('hltMu5L2Mu2JpsiTrackMassFiltered')"),
   Mu5_L2Mu3_Jpsi_MU = cms.string("!triggerObjectMatchesByCollection('hltL3MuonCandidates').empty() && "+
                                  " triggerObjectMatchesByCollection('hltL3MuonCandidates').at(0).hasFilterLabel('hltMu5L2Mu3JpsiTrackMassFiltered')"),
)

LowPtTriggerFlagsEfficienciesProbe = cms.PSet(
   ## Mu + Track
   Mu5_Track2_Jpsi_TK = cms.string("!triggerObjectMatchesByCollection('hltMuTrackJpsiEffCtfTrackCands').empty() && "+
                                   " triggerObjectMatchesByCollection('hltMuTrackJpsiEffCtfTrackCands').at(0).hasFilterLabel('hltMu5Track2JpsiTrackMassFiltered')"),
   Mu5_Track3p5_Jpsi_TK = cms.string("!triggerObjectMatchesByCollection('hltMuTrackJpsiEffCtfTrackCands').empty() && "+
                                   " triggerObjectMatchesByCollection('hltMuTrackJpsiEffCtfTrackCands').at(0).hasFilterLabel('hltMu5Track3p5JpsiTrackMassFiltered')"),
   Mu7_Track7_Jpsi_TK = cms.string("!triggerObjectMatchesByCollection('hltMuTrackJpsiCtfTrackCands').empty() && "+
                                   " triggerObjectMatchesByCollection('hltMuTrackJpsiCtfTrackCands').at(0).hasFilterLabel('hltMu7Track7JpsiTrackMassFiltered')"),
   ## Mu + L2Mu
   Mu5_L2Mu0_L2 = cms.string("!triggerObjectMatchesByCollection('hltL2MuonCandidates').empty() && "+
                             " triggerObjectMatchesByCollection('hltL2MuonCandidates').at(0).hasFilterLabel('hltDiMuonL2PreFiltered0')"),
   Mu5_L2Mu2_Jpsi_L2 = cms.string("!triggerObjectMatchesByCollection('hltL2MuonCandidates').empty() && "+
                             " triggerObjectMatchesByCollection('hltL2MuonCandidates').at(0).hasFilterLabel('hltMu5L2Mu2JpsiTrackMassFiltered')"),
   Mu5_L2Mu3_Jpsi_L2 = cms.string("!triggerObjectMatchesByCollection('hltL2MuonCandidates').empty() && "+
                             " triggerObjectMatchesByCollection('hltL2MuonCandidates').at(0).hasFilterLabel('hltMu5L2Mu3JpsiTrackMassFiltered')"),
)

LowPtTriggerFlagsEfficiencies = cms.PSet(LowPtTriggerFlagsEfficienciesTag,LowPtTriggerFlagsEfficienciesProbe)

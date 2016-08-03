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
    relTkIso  = cms.string("isolationR03.sumPt/pt"),
    ecalIso = cms.string("isolationR03.emEt"),
    relEcalIso = cms.string("isolationR03.emEt/pt"),
    hcalIso = cms.string("isolationR03.hadEt"),
    relHcalIso = cms.string("isolationR03.hadEt/pt"),
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
   # tkExpHitIn  = cms.string("? track.isNull ? 0 : track.trackerExpectedHitsInner().numberOfLostHits"),
  #  tkExpHitOut = cms.string("? track.isNull ? 0 : track.trackerExpectedHitsOuter().numberOfLostHits"),
    tkExpHitIn = cms.string("? track.isNull ? 0 : track.hitPattern.numberOfLostHits('MISSING_INNER_HITS')"), #reco::HitPattern::
    tkExpHitOut = cms.string("? track.isNull ? 0 : track.hitPattern.numberOfLostHits('MISSING_OUTER_HITS')"), #reco::HitPattern::
    #  tkHitFract  = cms.string("? track.isNull ? 0 : track.hitPattern.numberOfValidHits/(track.hitPattern.numberOfValidHits+track.hitPattern.numberOfLostHits+track.trackerExpectedHitsInner().numberOfLostHits+track.trackerExpectedHitsOuter().numberOfLostHits)"),
    tkHitFract  = cms.string("? track.isNull ? 0 : track.hitPattern.numberOfValidHits/(track.hitPattern.numberOfValidHits+track.hitPattern.numberOfLostHits('TRACK_HITS')+track.hitPattern.numberOfLostHits('MISSING_INNER_HITS')+ track.hitPattern.numberOfLostHits('MISSING_OUTER_HITS') )"),
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
    chi2LocMom = cms.string("combinedQuality.chi2LocalMomentum"),
    chi2LocPos = cms.string("combinedQuality.chi2LocalPosition"),
    glbTrackProb = cms.string("combinedQuality.glbTrackProbability"),
)
StaOnlyVariables = cms.PSet(
    staQoverP      = cms.string("? outerTrack.isNull() ? 0 : outerTrack.qoverp"),
    staQoverPerror = cms.string("? outerTrack.isNull() ? 0 : outerTrack.qoverpError"),
    staValidStations = cms.string("? outerTrack.isNull() ? -1 : outerTrack.hitPattern.muonStationsWithValidHits()"),
)
L1Variables = cms.PSet(
    l1pt = cms.string("? userCand('muonL1Info').isNull ? 0 : userCand('muonL1Info').pt"),
    l1eta  = cms.string("? userCand('muonL1Info').isNull ? 0 : userCand('muonL1Info').eta"),
    l1phi  = cms.string("? userCand('muonL1Info').isNull ? 0 : userCand('muonL1Info').phi"),
    l1q  = cms.string("userInt('muonL1Info:quality')"),
    l1dr = cms.string("userFloat('muonL1Info:deltaR')"),
    l1dphi = cms.string("userFloat('muonL1Info:deltaPhi')"),
    #l1ptByQ = cms.string("? userCand('muonL1Info:ByQ').isNull ? 0 : userCand('muonL1Info:ByQ').pt"),
    #l1qByQ  = cms.string("userInt('muonL1Info:qualityByQ')"),
    #l1drByQ = cms.string("userFloat('muonL1Info:deltaRByQ')"),
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
AllVariables = cms.PSet(KinematicVariables, IsolationVariables, MuonIDVariables, MuonCaloVariables, TrackQualityVariables, GlobalTrackQualityVariables, StaOnlyVariables, L1Variables, L2Variables, L3Variables)

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
    Loose       = cms.string("isLooseMuon()"),
    Medium      = cms.string("isPFMuon && innerTrack.validFraction >= 0.8 && ( isGlobalMuon && globalTrack.normalizedChi2 < 3 && combinedQuality.chi2LocalPosition < 12 && combinedQuality.trkKink < 20 && segmentCompatibility >= 0.303 || segmentCompatibility >= 0.451 )"),
    Medium2016      = cms.string("isPFMuon && innerTrack.validFraction >= 0.49 && ( isGlobalMuon && globalTrack.normalizedChi2 < 3 && combinedQuality.chi2LocalPosition < 12 && combinedQuality.trkKink < 20 && segmentCompatibility >= 0.303 || segmentCompatibility >= 0.451 )"),
    HighPt = cms.string("isGlobalMuon && isTrackerMuon &&  globalTrack.hitPattern.numberOfValidMuonHits > 0 && "+
                        "numberOfMatchedStations > 1 && track.hitPattern.numberOfValidPixelHits > 0 && "+
                        "track.hitPattern.trackerLayersWithMeasurement > 5 && abs(dB) < 0.2 && "+
                        "(tunePMuonBestTrack.ptError / tunePMuonBestTrack.pt) < 0.3"),
    HWWID =  cms.string("( ((isGlobalMuon() && "
                        "    globalTrack.normalizedChi2 <10 &&" +
                        "    globalTrack.hitPattern.numberOfValidMuonHits > 0 && " + 
                        "    numberOfMatches > 1 ) || " + 
                        "   (isTrackerMuon() && muonID('TMLastStationTight')) ) && " +
                        " isPFMuon && "+
                        " combinedQuality.trkKink < 20 &&" +
                        " innerTrack.hitPattern.trackerLayersWithMeasurement > 5 &&" +
                        " innerTrack.hitPattern.numberOfValidPixelHits > 0 && " + 
                        " abs(track.ptError / pt) < 0.10 )"),
    MuIDForOutsideInTk = cms.string("isStandAloneMuon && outerTrack.pt > 10 && outerTrack.hitPattern.muonStationsWithValidHits() >= 2"),
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
   #Mu40      = cms.string("!triggerObjectMatchesByPath('HLT_Mu40_v*',1,0).empty()"),
   #Mu40_eta2p1 = cms.string("!triggerObjectMatchesByPath('HLT_Mu40_eta2p1_v*',1,0).empty()"),
   #IsoMu15   = cms.string("!triggerObjectMatchesByPath('HLT_IsoMu15_v*',1,0).empty()"),
   #IsoMu15_eta2p1 = cms.string("!triggerObjectMatchesByPath('HLT_IsoMu15_eta2p1_v*',1,0).empty()"),
   #IsoMu17   = cms.string("!triggerObjectMatchesByPath('HLT_IsoMu17_v*',1,0).empty()"),
   #IsoMu24   = cms.string("!triggerObjectMatchesByPath('HLT_IsoMu24_v*',1,0).empty()"),
   #IsoMu24_eta2p1   = cms.string("!triggerObjectMatchesByPath('HLT_IsoMu24_eta2p1_v*',1,0).empty()"),
   #IsoMu30   = cms.string("!triggerObjectMatchesByPath('HLT_IsoMu30_v*',1,0).empty()"),
   #IsoMu30_eta2p1   = cms.string("!triggerObjectMatchesByPath('HLT_IsoMu30_eta2p1_v*',1,0).empty()"),
   #Mu15orMu17orMu20orMu24orMu30orMu40   = cms.string("!triggerObjectMatchesByPath('HLT_Mu15_v*',1,0).empty() || !triggerObjectMatchesByPath('HLT_IsoMu17_v*',1,0).empty() || !triggerObjectMatchesByPath('HLT_Mu20_v*',1,0).empty() || !triggerObjectMatchesByPath('HLT_IsoMu20_v*',1,0).empty() || !triggerObjectMatchesByPath('HLT_Mu24_v*',1,0).empty() || !triggerObjectMatchesByPath('HLT_Mu24_eta2p1_v*',1,0).empty() || !triggerObjectMatchesByPath('HLT_IsoMu24_eta2p1_v*',1,0).empty() || !triggerObjectMatchesByPath('HLT_IsoMu24_v*',1,0).empty() || !triggerObjectMatchesByPath('HLT_IsoMu30_eta2p1_v*',1,0).empty() || !triggerObjectMatchesByPath('HLT_IsoMu30_eta2p1_v*',1,0).empty() || !triggerObjectMatchesByPath('HLT_Mu30_v*',1,0).empty() || !triggerObjectMatchesByPath('HLT_Mu30_eta2p1_v*',1,0).empty() || !triggerObjectMatchesByPath('HLT_Mu40_v*',1,0).empty() || !triggerObjectMatchesByPath('HLT_Mu40_eta2p1_v*',1,0).empty()"),
   #DoubleMu5 = cms.string("!triggerObjectMatchesByPath('HLT_DoubleMu5_v*',1,0).empty()"),
   #DoubleMu7 = cms.string("!triggerObjectMatchesByPath('HLT_DoubleMu7_v*',1,0).empty()"),
   #DoubleMu13Mu8_Mu13 = cms.string("!triggerObjectMatchesByPath('HLT_Mu13_Mu8_v*',1,0).empty()"),
   #DoubleMu13Mu8_Mu8 = cms.string("!triggerObjectMatchesByFilter('hltDiMuonL3PreFiltered8').empty() || !triggerObjectMatchesByFilter('hltDiMuonL3p5PreFiltered8').empty()"),

   ## Heavily prescaled but still useful   
   #Mu17 = cms.string("!triggerObjectMatchesByPath('HLT_Mu17_v*',1,0).empty()"),
   #Mu8  = cms.string("!triggerObjectMatchesByPath('HLT_Mu8_v*',1,0).empty()"),

   # 2011 version
   # DoubleMu17Mu8_Mu17 = cms.string("!triggerObjectMatchesByPath('HLT_Mu17_Mu8_v*',1,0).empty()"),
   # DoubleMu17Mu8_Mu8 = cms.string("!triggerObjectMatchesByFilter('hltDiMuonL3PreFiltered8').empty() || !triggerObjectMatchesByFilter('hltDiMuonL3p5PreFiltered8').empty()"),
   # DoubleMu17TkMu8_TkMu8 = cms.string("!triggerObjectMatchesByPath('HLT_Mu17_TkMu8_v*',1,0).empty()"),
   # DoubleMu17TkMu8_Mu17 = cms.string("!triggerObjectMatchesByFilter('hltL3Mu17FromDiMuonFiltered17').empty()")

   # 2015 version
   Mu17             = cms.string("!triggerObjectMatchesByPath('HLT_Mu17_v*',1,0).empty()"),
   Mu17_IsoTrkVVL   = cms.string("!triggerObjectMatchesByPath('HLT_Mu17_TrkIsoVVL_v*',1,0).empty()"),
   Mu20             = cms.string("!triggerObjectMatchesByPath('HLT_Mu20_v*',1,0).empty()"),
   Mu50             = cms.string("!triggerObjectMatchesByPath('HLT_Mu50_v*',1,0).empty()"),
   Mu45_eta2p1      = cms.string("!triggerObjectMatchesByPath('HLT_Mu45_eta2p1_v*',1,0).empty()"),
   OldIsoMu18       = cms.string("!triggerObjectMatchesByPath('HLT_OldIsoMu18_v*',1,0).empty()"),
   IsoMu18          = cms.string("!triggerObjectMatchesByPath('HLT_IsoMu18_v*',1,0).empty()"),
   IsoMu20          = cms.string("!triggerObjectMatchesByPath('HLT_IsoMu20_v*',1,0).empty()"),
   IsoMu22          = cms.string("!triggerObjectMatchesByPath('HLT_IsoMu22_v*',1,0).empty()"),
   IsoMu24          = cms.string("!triggerObjectMatchesByPath('HLT_IsoMu24_v*',1,0).empty()"),
   IsoMu27          = cms.string("!triggerObjectMatchesByPath('HLT_IsoMu27_v*',1,0).empty()"),
   IsoMu17_eta2p1   = cms.string("!triggerObjectMatchesByPath('HLT_IsoMu17_eta2p1_v*',1,0).empty()"),
   IsoMu20_eta2p1   = cms.string("!triggerObjectMatchesByPath('HLT_IsoMu20_eta2p1_v*',1,0).empty()"),
   IsoMu24_eta2p1   = cms.string("!triggerObjectMatchesByPath('HLT_IsoMu24_eta2p1_v*',1,0).empty()"),
   OldIsoTkMu18     = cms.string("!triggerObjectMatchesByPath('HLT_OldIsoTkMu18_v*',1,0).empty()"),
   IsoTkMu18        = cms.string("!triggerObjectMatchesByPath('HLT_IsoTkMu18_v*',1,0).empty()"),
   IsoTkMu20        = cms.string("!triggerObjectMatchesByPath('HLT_IsoTkMu20_v*',1,0).empty()"),
   IsoTkMu22        = cms.string("!triggerObjectMatchesByPath('HLT_IsoTkMu22_v*',1,0).empty()"),
   IsoTkMu24        = cms.string("!triggerObjectMatchesByPath('HLT_IsoTkMu24_v*',1,0).empty()"),
   IsoTkMu27        = cms.string("!triggerObjectMatchesByPath('HLT_IsoTkMu27_v*',1,0).empty()"),
   IsoTkMu20_eta2p1 = cms.string("!triggerObjectMatchesByPath('HLT_IsoTkMu20_eta2p1_v*',1,0).empty()"),
   IsoTkMu24_eta2p1 = cms.string("!triggerObjectMatchesByPath('HLT_IsoTkMu24_eta2p1_v*',1,0).empty()"),
   HLT_TkMu50 = cms.string("!triggerObjectMatchesByPath('HLT_TkMu50_v*',1,0).empty()"),

   # To take care of the presence or not of the dz filter, we are requiring three flags for each of the main un-prescaled DoubleMuon Triggers
   # Mu17 leg, Mu8 leg, and (Mu17 leg && fired path)

   #DoubleMu17Mu8_Mu17 = cms.string("!triggerObjectMatchesByPath('HLT_Mu17_Mu8_v*',1,0).empty() && (!triggerObjectMatchesByFilter('hltL3fL1DoubleMu10MuOpenL1f0L2f10L3Filtered17').empty() || !triggerObjectMatchesByFilter('hltL3fL1DoubleMu10MuOpenOR3p5L1f0L2f10L3Filtered17').empty())"),
   ##DoubleMu17Mu8_Mu8 = cms.string("!triggerObjectMatchesByPath('HLT_Mu17_Mu8_v*',1,0).empty() && (!triggerObjectMatchesByFilter('hltL3pfL1DoubleMu10MuOpenL1f0L2pf0L3PreFiltered8').empty() || !triggerObjectMatchesByFilter('hltL3pfL1DoubleMu10MuOpenOR3p5L1f0L2pf0L3PreFiltered8').empty())"),
   #DoubleMu17Mu8_Mu17leg = cms.string("!triggerObjectMatchesByFilter('hltL3fL1DoubleMu10MuOpenL1f0L2f10L3Filtered17').empty() || !triggerObjectMatchesByFilter('hltL3fL1DoubleMu10MuOpenOR3p5L1f0L2f10L3Filtered17').empty()"),
   #DoubleMu17Mu8_Mu8leg = cms.string("!triggerObjectMatchesByFilter('hltL3pfL1DoubleMu10MuOpenL1f0L2pf0L3PreFiltered8').empty() || !triggerObjectMatchesByFilter('hltL3pfL1DoubleMu10MuOpenOR3p5L1f0L2pf0L3PreFiltered8').empty()"), 


   #DoubleMu17TkMu8_Mu17 = cms.string("!triggerObjectMatchesByPath('HLT_Mu17_TkMu8_v*',1,0).empty() && (!triggerObjectMatchesByFilter('hltL3fL1sMu10MuOpenL1f0L2f10L3Filtered17').empty() || !triggerObjectMatchesByFilter('hltL3fL1sMu10MuOpenOR3p5L1f0L2f10L3Filtered17').empty())"),
   #DoubleMu17TkMu8_TkMu8 = cms.string("!triggerObjectMatchesByPath('HLT_Mu17_TkMu8_v*',1,0).empty() && !triggerObjectMatchesByFilter('hltDiMuonGlbFiltered17TrkFiltered8').empty()"),
   #DoubleMu17TkMu8_Mu17leg = cms.string("!triggerObjectMatchesByFilter('hltL3fL1sMu10MuOpenL1f0L2f10L3Filtered17').empty() || !triggerObjectMatchesByFilter('hltL3fL1sMu10MuOpenOR3p5L1f0L2f10L3Filtered17').empty()"),
   #DoubleMu17TkMu8_TkMu8leg = cms.string("!triggerObjectMatchesByFilter('hltDiMuonGlbFiltered17TrkFiltered8').empty()"),

   #DoubleMu17TkMu8NoDZ_Mu17 = cms.string("!triggerObjectMatchesByPath('HLT_Mu17_TkMu8_NoDZ_v*',1,0).empty() && !triggerObjectMatchesByFilter('hltL3fL1sMu10MuOpenOR3p5L1f0L2f10L3Filtered17').empty()"),
   #DoubleMu13Mu8NoDZ_Mu13 = cms.string("!triggerObjectMatchesByPath('HLT_Mu13_Mu8_NoDZ_v*',1,0).empty() && !triggerObjectMatchesByFilter('hltL3fL1DoubleMu10MuOpenOR3p5L1f0L2f10L3Filtered13').empty()"),
   #DoubleMu13Mu8NoDZ_Mu8leg = cms.string("!triggerObjectMatchesByFilter('hltL3pfL1DoubleMu10MuOpenOR3p5L1f0L2pf0L3PreFiltered8').empty()"),

   # 2015 version
   DoubleIsoMu17Mu8dZ_Mu17leg = cms.string("!triggerObjectMatchesByPath('HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*',1,0).empty() && !triggerObjectMatchesByFilter('hltL3fL1sDoubleMu114L1f0L2f10OneMuL3Filtered17').empty()"),
   DoubleIsoMu17Mu8_IsoMu17leg = cms.string("!triggerObjectMatchesByFilter('hltDiMuonGlb17Glb8RelTrkIsoFiltered0p4').empty() && !triggerObjectMatchesByFilter('hltL3fL1sDoubleMu114L1f0L2f10OneMuL3Filtered17').empty()"),
   DoubleIsoMu17Mu8_Mu17leg = cms.string("!triggerObjectMatchesByFilter('hltL3fL1sDoubleMu114L1f0L2f10OneMuL3Filtered17').empty()"),
   DoubleIsoMu17Mu8_IsoMu8leg = cms.string("!triggerObjectMatchesByFilter('hltDiMuonGlb17Glb8RelTrkIsoFiltered0p4').empty() && (!triggerObjectMatchesByFilter('hltL3pfL1sDoubleMu114L1f0L2pf0L3PreFiltered8').empty()||!triggerObjectMatchesByFilter('hltL2pfL1sDoubleMu114ORDoubleMu125L1f0L2PreFiltered0').empty())"),
   DoubleIsoMu17Mu8_Mu8leg = cms.string("(!triggerObjectMatchesByFilter('hltL3pfL1sDoubleMu114L1f0L2pf0L3PreFiltered8').empty()||!triggerObjectMatchesByFilter('hltL2pfL1sDoubleMu114ORDoubleMu125L1f0L2PreFiltered0').empty())"),

   DoubleIsoMu17TkMu8dZ_Mu17 = cms.string("!triggerObjectMatchesByPath('HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*',1,0).empty() && !triggerObjectMatchesByFilter('hltL3fL1sDoubleMu114L1f0L2f10L3Filtered17').empty()"),
   DoubleIsoMu17TkMu8_IsoMu17leg = cms.string("!triggerObjectMatchesByFilter('hltDiMuonGlb17Trk8RelTrkIsoFiltered0p4').empty() && !triggerObjectMatchesByFilter('hltL3fL1sDoubleMu114L1f0L2f10L3Filtered17').empty()"),
   DoubleIsoMu17TkMu8_Mu17leg = cms.string("!triggerObjectMatchesByFilter('hltL3fL1sDoubleMu114L1f0L2f10L3Filtered17').empty()"),
   DoubleIsoMu17TkMu8_IsoMu8leg = cms.string("!triggerObjectMatchesByFilter('hltDiMuonGlb17Trk8RelTrkIsoFiltered0p4').empty() && !triggerObjectMatchesByFilter('hltDiMuonGlbFiltered17TrkFiltered8').empty()"),
   DoubleIsoMu17TkMu8_TkMu8leg = cms.string("!triggerObjectMatchesByFilter('hltDiMuonGlbFiltered17TrkFiltered8').empty()"),
     
                              
   #monitoring of the Mu30TkMu11 path
   DoubleMu30TkMu11 = cms.string("!triggerObjectMatchesByPath('HLT_Mu30_TkMu11_v*',1,0).empty()"),
   DoubleMu30TkMu11_Mu30leg = cms.string("!triggerObjectMatchesByFilter('hltL3fL1sMu16orMu25L1f0L2f25L3Filtered30').empty()"),
   DoubleMu30TkMu11_TkMu11leg = cms.string("!triggerObjectMatchesByFilter('hltDiMuonGlbFiltered30TrkFiltered11').empty()"),

)
HighPtTriggerFlagsDebug = cms.PSet(
   # --- the ones commented out don't save tags ---
   #L1DoubleMu10MuOpenL1Filtered0 = cms.string("!triggerObjectMatchesByFilter('hltL1DoubleMu10MuOpenL1Filtered0').empty()"),
   #L1DoubleMu10MuOpenOR3p5L1Filtered0 = cms.string("!triggerObjectMatchesByFilter('hltL1DoubleMu10MuOpenOR3p5L1Filtered0').empty()"),
   #L1fL1sDoubleMu10MuOpenL1Filtered0 = cms.string("!triggerObjectMatchesByFilter('hltL1fL1sDoubleMu10MuOpenL1Filtered0').empty()"),
   #L1fL1sDoubleMu10MuOpenOR3p5L1Filtered0  = cms.string("!triggerObjectMatchesByFilter('hltL1fL1sDoubleMu10MuOpenOR3p5L1Filtered0').empty()"),
   #L1fL1sMu16Eta2p1L1Filtered0 = cms.string("!triggerObjectMatchesByFilter('hltL1fL1sMu16Eta2p1L1Filtered0').empty()"),
   #L1fL1sMu16L1Filtered0 = cms.string("!triggerObjectMatchesByFilter('hltL1fL1sMu16L1Filtered0').empty()"),
   L1sMu16 = cms.string("!triggerObjectMatchesByFilter('hltL1sSingleMu16').empty()"),
   L1sMu16Eta2p1 = cms.string("!triggerObjectMatchesByFilter('hltL1sL1SingleMu16eta2p1').empty()"),
   L1sMu18 = cms.string("!triggerObjectMatchesByFilter('hltL1sSingleMu18').empty()"),
   L1sMu20 = cms.string("!triggerObjectMatchesByFilter('hltL1sSingleMu20').empty()"),
   L1sMu25Eta2p1 = cms.string("!triggerObjectMatchesByFilter('hltL1sL1SingleMu25eta2p1').empty()"),
   L2pfL1sDoubleMu114L1f0L2PreFiltered0 = cms.string("!triggerObjectMatchesByFilter('hltL2pfL1sDoubleMu114L1f0L2PreFiltered0').empty()"),
   L2fL1sDoubleMu114L1f0L2Filtered10OneMu = cms.string("!triggerObjectMatchesByFilter('hltL2fL1sDoubleMu114L1f0L2Filtered10OneMu').empty()"),
   L3pfL1sDoubleMu114L1f0L2pf0L3PreFiltered8 = cms.string("!triggerObjectMatchesByFilter('hltL3pfL1sDoubleMu114L1f0L2pf0L3PreFiltered8').empty()"),
   L3fL1sDoubleMu114L1f0L2f10OneMuL3Filtered17 = cms.string("!triggerObjectMatchesByFilter('hltL3fL1sDoubleMu114L1f0L2f10OneMuL3Filtered17').empty()"),
   DiMuonGlb17Glb8RelTrkIsoFiltered0p4 = cms.string("!triggerObjectMatchesByFilter('hltDiMuonGlb17Glb8RelTrkIsoFiltered0p4').empty()"),
   L2fL1sDoubleMu114L1f0OneMuL2Filtered10 = cms.string("!triggerObjectMatchesByFilter('hltL2fL1sDoubleMu114L1f0OneMuL2Filtered10').empty()"),
   L3fL1sDoubleMu114L1f0L2f10L3Filtered17 = cms.string("!triggerObjectMatchesByFilter('hltL3fL1sDoubleMu114L1f0L2f10L3Filtered17').empty()"),
   DiMuonGlbFiltered17TrkFiltered8 = cms.string("!triggerObjectMatchesByFilter('hltDiMuonGlbFiltered17TrkFiltered8').empty()"),
   DiMuonGlb17Trk8RelTrkIsoFiltered0p4 = cms.string("!triggerObjectMatchesByFilter('L2pfL1DoubleMu10MuOpenOR3p5L1f0L2PreFiltered0').empty()"),
   L3fL1sMu16Eta2p1L1f0L2f16QL3Filtered24Q = cms.string("!triggerObjectMatchesByFilter('hltL3fL1sMu16Eta2p1L1f0L2f16QL3Filtered24Q').empty()"),
   L3fL1sMu16L1f0L2f16QL3Filtered24Q = cms.string("!triggerObjectMatchesByFilter('hltL3fL1sMu16L1f0L2f16QL3Filtered24Q').empty()"),
   L3fL1sMu16L1f0L2f10QL3Filtered18Q = cms.string("!triggerObjectMatchesByFilter('hltL3fL1sMu16L1f0L2f10QL3Filtered18Q').empty()"),
   L3fL1sMu16L1f0L2f10QL3Filtered18QL3pfecalIsoRhoFilteredEB0p11EE0p08 = cms.string("!triggerObjectMatchesByFilter('hltL3fL1sMu16L1f0L2f10QL3Filtered18QL3pfecalIsoRhoFilteredEB0p11EE0p08').empty()"),
   L3fL1sMu16L1f0L2f10QL3Filtered18QL3pfhcalIsoRhoFilteredHB0p21HE0p22 = cms.string("!triggerObjectMatchesByFilter('hltL3fL1sMu16L1f0L2f10QL3Filtered18QL3pfhcalIsoRhoFilteredHB0p21HE0p22').empty()"),
   L3crIsoL1sMu16L1f0L2f10QL3f18QL3trkIsoFiltered0p09 = cms.string("!triggerObjectMatchesByFilter('hltL3crIsoL1sMu16L1f0L2f10QL3f18QL3trkIsoFiltered0p09').empty()"),
   L3fL1sMu16L1f0L2f10QL3Filtered18QL3pfecalOldIsoRhoFilteredEB0p11EE0p08 = cms.string("!triggerObjectMatchesByFilter('hltL3fL1sMu16L1f0L2f10QL3Filtered18QL3pfecalOldIsoRhoFilteredEB0p11EE0p08').empty()"),
   L3fL1sMu16L1f0L2f10QL3Filtered18QL3pfhcalOldIsoRhoFilteredHB0p21HE0p22 = cms.string("!triggerObjectMatchesByFilter('hltL3fL1sMu16L1f0L2f10QL3Filtered18QL3pfhcalOldIsoRhoFilteredHB0p21HE0p22').empty()"),
   L3crIsoL1sMu16L1f0L2f10QL3f18QL3OldCaloIsotrkIsoFiltered0p09 = cms.string("!triggerObjectMatchesByFilter('hltL3crIsoL1sMu16L1f0L2f10QL3f18QL3OldCaloIsotrkIsoFiltered0p09').empty()"),
   L3fL1sMu16f0TkFiltered18Q = cms.string("!triggerObjectMatchesByFilter('hltL3fL1sMu16f0TkFiltered18Q').empty()"),
   #add the control path HLT_Mu17_TrkIsoVVL
   hltL2fL1sMu10lqL1f0L2Filtered10 = cms.string("!triggerObjectMatchesByFilter('hltL2fL1sMu10lqL1f0L2Filtered10').empty()"),
   hltL3fL1sMu10lqL1f0L2f10L3Filtered17 = cms.string("!triggerObjectMatchesByFilter('hltL3fL1sMu10lqL1f0L2f10L3Filtered17').empty()"),
   hltL3fL1sMu1lqL1f0L2f10L3Filtered17TkIsoFiltered0p4 = cms.string("!triggerObjectMatchesByFilter('hltL3fL1sMu1lqL1f0L2f10L3Filtered17TkIsoFiltered0p4').empty()"),

   #add the muon 20 isolation variables
   L2fL1sMu18L1f0L2Filtered10Q = cms.string("!triggerObjectMatchesByFilter('hltL2fL1sMu18L1f0L2Filtered10Q').empty()"),
   hltL3fL1sMu18L1f0L2f10QL3Filtered20Q = cms.string("!triggerObjectMatchesByFilter('hltL3fL1sMu18L1f0L2f10QL3Filtered20Q').empty()"),
   hltL3fL1sMu18L1f0L2f10QL3Filtered20QL3pfecalIsoRhoFilteredEB0p11EE0p08 = cms.string("!triggerObjectMatchesByFilter('hltL3fL1sMu18L1f0L2f10QL3Filtered20QL3pfecalIsoRhoFilteredEB0p11EE0p08').empty()"),
   hltL3fL1sMu18L1f0L2f10QL3Filtered20QL3pfhcalIsoRhoFilteredHB0p21HE0p22 = cms.string("!triggerObjectMatchesByFilter('hltL3fL1sMu18L1f0L2f10QL3Filtered20QL3pfhcalIsoRhoFilteredHB0p21HE0p22').empty()"),
   hltL3crIsoL1sMu18L1f0L2f10QL3f20QL3trkIsoFiltered0p09 = cms.string("!triggerObjectMatchesByFilter('hltL3crIsoL1sMu18L1f0L2f10QL3f20QL3trkIsoFiltered0p09').empty()"),
   hltL3fL1sMu18f0TkFiltered20Q = cms.string("!triggerObjectMatchesByFilter('hltL3fL1sMu18f0TkFiltered20Q').empty()"),
   hltL3fL1sMu18f0TkFiltered20QL3pfecalIsoRhoFilteredEB0p11EE0p08 = cms.string("!triggerObjectMatchesByFilter('hltL3fL1sMu18f0TkFiltered20QL3pfecalIsoRhoFilteredEB0p11EE0p08').empty()"),
   hltL3fL1sMu18f0TkFiltered20QL3pfhcalIsoRhoFilteredHB0p21HE0p22 = cms.string("!triggerObjectMatchesByFilter('hltL3fL1sMu18f0TkFiltered20QL3pfhcalIsoRhoFilteredHB0p21HE0p22').empty()"),
   hltL3fL1sMu18L1f0Tkf20QL3trkIsoFiltered0p09 = cms.string("!triggerObjectMatchesByFilter('hltL3fL1sMu18L1f0Tkf20QL3trkIsoFiltered0p09').empty()"),

   #matching to muon leg of HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v*
   Mu23_TrkIsoVVL = cms.string("!triggerObjectMatchesByFilter('hltMu23TrkIsoVVLEle8CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered23').empty()"),
   
)


LowPtTriggerFlagsPhysics = cms.PSet(
    ## L2 filters 
    Dimuon16_L1L2 = cms.string(         "!triggerObjectMatchesByCollection('hltL2MuonCandidates').empty() && "+ 
                                        "(triggerObjectMatchesByCollection('hltL2MuonCandidates').at(0).hasFilterLabel('hltL2fL1sL1DoubleMu10MuOpenL1f0L2PreFiltered0') || "+
                                        " triggerObjectMatchesByCollection('hltL2MuonCandidates').at(0).hasFilterLabel('hltL2fL1sL1DoubleMu100dEtaMax1p8IorDoubleMu114L1f0L2PreFiltered0'))"
                                        ),
    Dimuon10_L1L2 = cms.string(         "!triggerObjectMatchesByCollection('hltL2MuonCandidates').empty() && "+
                                        "(triggerObjectMatchesByCollection('hltL2MuonCandidates').at(0).hasFilterLabel('hltL2fL1sL1DoubleMu0er16NoOSL1f0L2PreFiltered0') || "+
                                        " triggerObjectMatchesByCollection('hltL2MuonCandidates').at(0).hasFilterLabel('hltL2fL1sL1DoubleMu0er16IorDoubleMu0er16OSL1f0L2PreFiltered0'))"
                                        ),
    ## DoubleMu triggers filters 
    DoubleMu4_L1L2 = cms.string(        "!triggerObjectMatchesByCollection('hltL2MuonCandidates').empty() && "+
                                        "(triggerObjectMatchesByCollection('hltL2MuonCandidates').at(0).hasFilterLabel('hltL2fL1sL1DoubleMu0er16ORL1DoubleMu10MuOpenL1f0L2PreFiltered0') || "+
                                        " triggerObjectMatchesByCollection('hltL2MuonCandidates').at(0).hasFilterLabel('hltL2fL1sL1DoubleMu0er16IorDoubleMu0er16OSIorL1DoubleMu10MuOpenL1f0L2PreFiltered0'))"
                                        ),

    ## L3 Mu 
    Mu_L3 = cms.string(                 "!triggerObjectMatchesByCollection('hltL3MuonCandidates').empty()"
                                        ),
    ## Vertexing filters 
    Dimuon6_Jpsi_NoVertexing = cms.string("!triggerObjectMatchesByCollection('hltL3MuonCandidates').empty() && "+
                                          " triggerObjectMatchesByCollection('hltL3MuonCandidates').at(0).hasFilterLabel('hltDimuon6JpsiL3Filtered')"),
    Dimuon0er16_Jpsi_NoVertexing = cms.string("!triggerObjectMatchesByCollection('hltL3MuonCandidates').empty() && "+ 
                                              " triggerObjectMatchesByCollection('hltL3MuonCandidates').at(0).hasFilterLabel('hltDimuon0JpsiOSL3Filtered')"),
    Dimuon0er16_Jpsi_NoOS_NoVertexing = cms.string("!triggerObjectMatchesByCollection('hltL3MuonCandidates').empty() && "+
                                                   " triggerObjectMatchesByCollection('hltL3MuonCandidates').at(0).hasFilterLabel('hltDimuon0JpsiNoOSL3Filtered')"),
    ## Final filters  
    Dimuon16_Jpsi = cms.string("!triggerObjectMatchesByCollection('hltL3MuonCandidates').empty() && "+
                               " triggerObjectMatchesByCollection('hltL3MuonCandidates').at(0).hasFilterLabel('hltDisplacedmumuFilterDimuon16Jpsi')"),
    Dimuon20_Jpsi = cms.string("!triggerObjectMatchesByCollection('hltL3MuonCandidates').empty() && "+
                               " triggerObjectMatchesByCollection('hltL3MuonCandidates').at(0).hasFilterLabel('hltDisplacedmumuFilterDimuon20Jpsi')"),
    Dimuon10_Jpsi_Barrel = cms.string("!triggerObjectMatchesByCollection('hltL3MuonCandidates').empty() && "+
                                      " triggerObjectMatchesByCollection('hltL3MuonCandidates').at(0).hasFilterLabel('hltDisplacedmumuFilterDimuon10JpsiBarrel')"),
)

LowPtTriggerFlagsEfficienciesTag = cms.PSet(
   ## Mu + Track
   ## 2012
   Mu5_Track2_Jpsi_MU = cms.string("!triggerObjectMatchesByCollection('hltL3MuonCandidates').empty() && "+
                                   " triggerObjectMatchesByCollection('hltL3MuonCandidates').at(0).hasFilterLabel('hltMu5Track2JpsiTrackMassFiltered')"),
   Mu5_Track3p5_Jpsi_MU = cms.string("!triggerObjectMatchesByCollection('hltL3MuonCandidates').empty() && "+
                                   " triggerObjectMatchesByCollection('hltL3MuonCandidates').at(0).hasFilterLabel('hltMu5Track3p5JpsiTrackMassFiltered')"),
   Mu7_Track7_Jpsi_MU = cms.string("!triggerObjectMatchesByCollection('hltL3MuonCandidates').empty() && "+
                                   " triggerObjectMatchesByCollection('hltL3MuonCandidates').at(0).hasFilterLabel('hltMu7Track7JpsiTrackMassFiltered')"),

   ## 2015
   Mu7p5_MU = cms.string("!triggerObjectMatchesByCollection('hltL3MuonCandidates').empty() && " + 
                         " triggerObjectMatchesByCollection('hltL3MuonCandidates').at(0).hasFilterLabel('hltL3fLMu7p5TrackL3Filtered7p5')"
                         ),

   Mu7p5_Track2_Jpsi_MU = cms.string("!triggerObjectMatchesByCollection('hltL3MuonCandidates').empty() && "+
                                     " triggerObjectMatchesByCollection('hltL3MuonCandidates').at(0).hasFilterLabel('hltMu7p5Track2JpsiTrackMassFiltered')"),
   Mu7p5_Track3p5_Jpsi_MU = cms.string("!triggerObjectMatchesByCollection('hltL3MuonCandidates').empty() && "+
                                       " triggerObjectMatchesByCollection('hltL3MuonCandidates').at(0).hasFilterLabel('hltMu7p5Track3p5JpsiTrackMassFiltered')"),
   Mu7p5_Track7_Jpsi_MU = cms.string("!triggerObjectMatchesByCollection('hltL3MuonCandidates').empty() && "+
                                     " triggerObjectMatchesByCollection('hltL3MuonCandidates').at(0).hasFilterLabel('hltMu7p5Track7JpsiTrackMassFiltered')"),

   ## Mu + L2Mu
   ## 2012
   Mu5_L2Mu0_MU = cms.string("!triggerObjectMatchesByCollection('hltL3MuonCandidates').empty() && "+
                             " triggerObjectMatchesByCollection('hltL3MuonCandidates').at(0).hasFilterLabel('hltMu5L2Mu0L3Filtered5')"),
   Mu5_L2Mu2_Jpsi_MU = cms.string("!triggerObjectMatchesByCollection('hltL3MuonCandidates').empty() && "+
                                  " triggerObjectMatchesByCollection('hltL3MuonCandidates').at(0).hasFilterLabel('hltMu5L2Mu2JpsiTrackMassFiltered')"),
   Mu5_L2Mu3_Jpsi_MU = cms.string("!triggerObjectMatchesByCollection('hltL3MuonCandidates').empty() && "+
                                  " triggerObjectMatchesByCollection('hltL3MuonCandidates').at(0).hasFilterLabel('hltMu5L2Mu3JpsiTrackMassFiltered')"),
   ## 2015
   Mu7p5_L2Mu2_Jpsi_MU = cms.string("!triggerObjectMatchesByCollection('hltL3MuonCandidates').empty() && "+
                                 " triggerObjectMatchesByCollection('hltL3MuonCandidates').at(0).hasFilterLabel('hltMu7p5L2Mu2JpsiTrackMassFiltered')"),

   ## Mu + TkMu
   ## 2015
   Mu25TkMu0Onia_L3_MU = cms.string(   "!triggerObjectMatchesByCollection('hltL3MuonCandidates').empty() && " +
                                       " triggerObjectMatchesByCollection('hltL3MuonCandidates').at(0).hasFilterLabel('hltL3fL1sMu16orMu20erorMu25L1f0L2f0L3Filtered25')"
                                       ),
   Mu16TkMu0Onia_L3_MU = cms.string(   "!triggerObjectMatchesByCollection('hltL3MuonCandidates').empty() && " + 
                                       " triggerObjectMatchesByCollection('hltL3MuonCandidates').at(0).hasFilterLabel('hltL3fL1sMu16orMu20erorMu16L1f0L2f0L3Filtered16')"
                                       ),

)

LowPtTriggerFlagsEfficienciesProbe = cms.PSet(

   ## Mu + Track
   ## 2012
   Mu5_Track2_Jpsi_TK = cms.string("!triggerObjectMatchesByCollection('hltMuTrackJpsiEffCtfTrackCands').empty() && "+
                                   " triggerObjectMatchesByCollection('hltMuTrackJpsiEffCtfTrackCands').at(0).hasFilterLabel('hltMu5Track2JpsiTrackMassFiltered')"),
   Mu5_Track3p5_Jpsi_TK = cms.string("!triggerObjectMatchesByCollection('hltMuTrackJpsiEffCtfTrackCands').empty() && "+
                                     " triggerObjectMatchesByCollection('hltMuTrackJpsiEffCtfTrackCands').at(0).hasFilterLabel('hltMu5Track3p5JpsiTrackMassFiltered')"),
   Mu7_Track7_Jpsi_TK = cms.string("!triggerObjectMatchesByCollection('hltMuTrackJpsiCtfTrackCands').empty() && "+
                                   " triggerObjectMatchesByCollection('hltMuTrackJpsiCtfTrackCands').at(0).hasFilterLabel('hltMu7Track7JpsiTrackMassFiltered')"),
    
   ## 2015
   Mu7p5_Track2_Jpsi_TK = cms.string("!triggerObjectMatchesByCollection('hltTracksIter').empty() && " + 
                                     " triggerObjectMatchesByCollection('hltTracksIter').at(0).hasFilterLabel('hltMu7p5Track2JpsiTrackMassFiltered')"),
   Mu7p5_Track3p5_Jpsi_TK = cms.string("!triggerObjectMatchesByCollection('hltTracksIter').empty() && " + 
                                       " triggerObjectMatchesByCollection('hltTracksIter').at(0).hasFilterLabel('hltMu7p5Track3p5JpsiTrackMassFiltered')"),
   Mu7p5_Track7_Jpsi_TK = cms.string("!triggerObjectMatchesByCollection('hltTracksIter').empty() && " + 
                                     " triggerObjectMatchesByCollection('hltTracksIter').at(0).hasFilterLabel('hltMu7p5Track7JpsiTrackMassFiltered')"),

   ## Mu + L2Mu
   ## 2012
   Mu5_L2Mu0_L2 = cms.string("!triggerObjectMatchesByCollection('hltL2MuonCandidates').empty() && "+
                             " triggerObjectMatchesByCollection('hltL2MuonCandidates').at(0).hasFilterLabel('hltDiMuonL2PreFiltered0')"),
   Mu5_L2Mu2_Jpsi_L2 = cms.string("!triggerObjectMatchesByCollection('hltL2MuonCandidates').empty() && "+
                                  " triggerObjectMatchesByCollection('hltL2MuonCandidates').at(0).hasFilterLabel('hltMu5L2Mu2JpsiTrackMassFiltered')"),
   Mu5_L2Mu3_Jpsi_L2 = cms.string("!triggerObjectMatchesByCollection('hltL2MuonCandidates').empty() && "+
                                  " triggerObjectMatchesByCollection('hltL2MuonCandidates').at(0).hasFilterLabel('hltMu5L2Mu3JpsiTrackMassFiltered')"),
   ## 2015
   Mu7p5_L2Mu2_L2 = cms.string("!triggerObjectMatchesByCollection('hltL2MuonCandidates').empty() && "+
                               " triggerObjectMatchesByCollection('hltL2MuonCandidates').at(0).hasFilterLabel('hltL2fDoubleMu2L2PreFiltered2')"
                                        ),
   Mu7p5_L2Mu2_Jpsi_L2 = cms.string("!triggerObjectMatchesByCollection('hltL2MuonCandidates').empty() && " + 
                                    " triggerObjectMatchesByCollection('hltL2MuonCandidates').at(0).hasFilterLabel('hltMu7p5L2Mu2JpsiTrackMassFiltered')"),

   ## Mu + TkMu 
   Mu25TkMu0Onia_L2 = cms.string(      "!triggerObjectMatchesByCollection('hltL2MuonCandidates').empty() && "+
                                       " triggerObjectMatchesByCollection('hltL2MuonCandidates').at(0).hasFilterLabel('hltL2fL1sMu16orMu20erorMu25L1f0L2Filtered0')"
                                       ),
   Mu25TkMu0Onia_L3 = cms.string(      "!triggerObjectMatchesByCollection('hltL3MuonCandidates').empty() && "+
                                       " triggerObjectMatchesByCollection('hltL3MuonCandidates').at(0).hasFilterLabel('hltL3fL1sMu16orMu20erorMu25L1f0L2f0L3Filtered25')"
                                       ),
   Mu25TkMu0Onia_TM = cms.string(      "!triggerObjectMatchesByCollection('hltGlbTrkMuonCands').empty()"
                                       ),
   Mu16TkMu0Onia_L2 = cms.string(      "!triggerObjectMatchesByCollection('hltL2MuonCandidates').empty() && "+
                                       " triggerObjectMatchesByCollection('hltL2MuonCandidates').at(0).hasFilterLabel('hltL2fL1sMu16orMu20erorMu16L1f0L2Filtered0')"
                                       ),
   Mu16TkMu0Onia_L3 = cms.string(      "!triggerObjectMatchesByCollection('hltL3MuonCandidates').empty() && "+
                                       " triggerObjectMatchesByCollection('hltL3MuonCandidates').at(0).hasFilterLabel('hltL3fL1sMu16orMu20erorMu16L1f0L2f0L3Filtered16')"
                                        ),

)

LowPtTriggerFlagsEfficiencies = cms.PSet(LowPtTriggerFlagsEfficienciesTag,LowPtTriggerFlagsEfficienciesProbe)

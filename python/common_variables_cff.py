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
    combRelIso = cms.string("(isolationR03.emEt + isolationR03.hadEt + isolationR03.sumPt)/pt"),
)
MuonIDVariables = cms.PSet(
    caloCompatibility = cms.string("? isCaloCompatibilityValid ? caloCompatibility : -1"),
    numberOfMatches   = cms.string("? isMatchesValid ? numberOfMatches : -1"),
    numberOfMatchedStations = cms.string("? isMatchesValid ? numberOfMatchedStations : -1"),
)
TrackQualityVariables = cms.PSet(
    dB          = cms.string("dB"),
    tkValidHits = cms.string("? track.isNull ? 0 : track.numberOfValidHits"),
    tkPixelHits = cms.string("? track.isNull ? 0 : track.hitPattern.numberOfValidPixelHits"),
    tkPixelLay  = cms.string("? track.isNull ? 0 : track.hitPattern.pixelLayersWithMeasurement"),
    tkExpHitIn  = cms.string("? track.isNull ? 0 : track.trackerExpectedHitsInner.numberOfLostHits"),
    tkExpHitInNew  =cms.string("userInt('expectedHitsMu:in')"),
    tkExpHitOutNew =cms.string("userInt('expectedHitsMu:out')")
    
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
AllVariables = cms.PSet(KinematicVariables, IsolationVariables, MuonIDVariables, TrackQualityVariables, L1Variables, L2Variables, L3Variables)

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
    TMOST = cms.string("muonID('TMOneStationTight')"),
     
    VBTF   = cms.string("numberOfMatchedStations > 1 && muonID('GlobalMuonPromptTight') && abs(dB) < 0.2 && "+
                        "track.numberOfValidHits > 10 && track.hitPattern.numberOfValidPixelHits > 0"),
    VBTFold   = cms.string("numberOfMatches > 1 && muonID('GlobalMuonPromptTight') && abs(dB) < 0.2 && "+
                        "track.numberOfValidHits > 10 && track.hitPattern.numberOfValidPixelHits > 0"),
    BMuons = cms.string("(isGlobalMuon && globalTrack.chi2 / globalTrack.ndof < 20. && "+ 
                        " globalTrack.hitPattern.numberOfValidMuonHits > 0 && isTrackerMuon && muonID('TrackerMuonArbitrated') )"+
                        "|| (isTrackerMuon && muonID('TrackerMuonArbitrated') )"),
    BMuQual = cms.string("((isGlobalMuon && globalTrack.chi2 / globalTrack.ndof < 20. && " + 
                        " globalTrack.hitPattern.numberOfValidMuonHits > 0 && isTrackerMuon && muonID('TrackerMuonArbitrated') )"+
                        "|| (isTrackerMuon && muonID('TrackerMuonArbitrated') ) ) && " +
                        "(track.numberOfValidHits > 10 && track.chi2 / track.ndof < 1.8 && muonID('TMOneStationTight') && "+
                        " track.hitPattern.pixelLayersWithMeasurement > 1 && abs(dB) < 3.0 && abs(track.dz) < 15.0)"),
    BBbarCSWithIP = cms.string("isGlobalMuon && isTrackerMuon "
                               +"&& track.numberOfValidHits > 10"
                               +"&& track.hitPattern.numberOfValidPixelHits > 1 "
                               +"&& track.normalizedChi2 < 2 && globalTrack.normalizedChi2 < 10 "
                               +"&& globalTrack.hitPattern.numberOfValidMuonHits > 0 && numberOfMatches > 1"
                               +"&& muonID('TrackerMuonArbitrated') && muonID('TMOneStationTight')"),


)

HighPtTriggerFlags = cms.PSet(
   Mu9       = cms.string("!triggerObjectMatchesByPath('HLT_Mu9').empty()"),
   Mu11      = cms.string("!triggerObjectMatchesByPath('HLT_Mu11').empty()"),
   Mu15      = cms.string("!triggerObjectMatchesByPath('HLT_Mu15_v1').empty()"),
   IsoMu9    = cms.string("!triggerObjectMatchesByPath('HLT_IsoMu9').empty()"),
   IsoMu13   = cms.string("!triggerObjectMatchesByPath('HLT_IsoMu13_v1').empty() || !triggerObjectMatchesByPath('HLT_IsoMu13_v2').empty() || !triggerObjectMatchesByPath('HLT_IsoMu13_v3').empty() || !triggerObjectMatchesByPath('HLT_IsoMu13_v4').empty()"),
   DoubleMu3 = cms.string("(!triggerObjectMatchesByPath('HLT_DoubleMu3').empty() || !triggerObjectMatchesByPath('HLT_DoubleMu3_v1').empty() || !triggerObjectMatchesByPath('HLT_DoubleMu3_v2').empty() )"),
   DoubleMu5 = cms.string("(!triggerObjectMatchesByPath('HLT_DoubleMu5').empty() || !triggerObjectMatchesByPath('HLT_DoubleMu5_v1').empty() || !triggerObjectMatchesByPath('HLT_DoubleMu5_v2').empty() )"),
)

LowPtTriggerFlagsPhysics = cms.PSet(
   #L1DoubleMuOpen       = cms.string("!triggerObjectMatchesByFilter('hltDoubleMuLevel1PathL1OpenFiltered').empty()"),
   #L1DoubleMuOpen_Tight = cms.string("!triggerObjectMatchesByFilter('hltL1DoubleMuOpenTightL1Filtered').empty()"),
   #L2DoubleMu0          = cms.string("!triggerObjectMatchesByFilter('hltDiMuonL2PreFiltered0').empty()"),
   DoubleMu0             = cms.string("!triggerObjectMatchesByFilter('hltDiMuonL3PreFiltered0').empty()"),
   DoubleMu0_Quarkonium  = cms.string("!triggerObjectMatchesByPath('HLT_DoubleMu0').empty() || "+
                                      "!triggerObjectMatchesByPath('HLT_DoubleMu0_Quarkonium').empty() || "+
                                      "!triggerObjectMatchesByPath('HLT_DoubleMu0_Quarkonium_v1').empty()"),

   ## Mu + Tracker Mu

   Mu0_TkMu0_Jpsi_MU = cms.string("!triggerObjectMatchesByCollection('hltL3MuonCandidates').empty() && "+
                                  " triggerObjectMatchesByCollection('hltL3MuonCandidates').at(0).hasFilterLabel('hltMu0TkMuJpsiTkMuMassFiltered')"),
   Mu0_TkMu0_Jpsi_TM = cms.string("!triggerObjectMatchesByCollection('hltMuTkMuJpsiTrackerMuonCands').empty() && "+
                                  " triggerObjectMatchesByCollection('hltMuTkMuJpsiTrackerMuonCands').at(0).hasFilterLabel('hltMu0TkMuJpsiTkMuMassFiltered')"),

   Mu0_TkMu0_OST_Jpsi_MU = cms.string("!triggerObjectMatchesByCollection('hltL3MuonCandidates').empty() && "+
                                      " triggerObjectMatchesByCollection('hltL3MuonCandidates').at(0).hasPathName('HLT_Mu0_TkMu0_OST_Jpsi')"),
   Mu0_TkMu0_OST_Jpsi_TK = cms.string("!triggerObjectMatchesByCollection('hltMuTkMuJpsiTrackerMuonCands').empty() && "+
                                      " triggerObjectMatchesByCollection('hltMuTkMuJpsiTrackerMuonCands').at(0).hasPathName('HLT_Mu0_TkMu0_OST_Jpsi')"),
  
   Mu0_TkMu0_OST_Tight_Jpsi_MU = cms.string("!triggerObjectMatchesByCollection('hltL3MuonCandidates').empty() && "+
                                            " (triggerObjectMatchesByCollection('hltL3MuonCandidates').at(0).hasPathName('HLT_Mu0_TkMu0_OST_Jpsi_Tight_v1') || "+
                                            "  triggerObjectMatchesByCollection('hltL3MuonCandidates').at(0).hasPathName('HLT_Mu0_TkMu0_OST_Jpsi_Tight_v2') || "+
                                            "  triggerObjectMatchesByCollection('hltL3MuonCandidates').at(0).hasPathName('HLT_Mu0_TkMu0_OST_Jpsi_Tight_v3') )"),
   Mu0_TkMu0_OST_Tight_Jpsi_TK = cms.string("!triggerObjectMatchesByCollection('hltMuTkMuJpsiTrackerMuonCands').empty() && "+
                                            " (triggerObjectMatchesByCollection('hltMuTkMuJpsiTrackerMuonCands').at(0).hasPathName('HLT_Mu0_TkMu0_OST_Jpsi_Tight_v1') || "+
                                            "  triggerObjectMatchesByCollection('hltMuTkMuJpsiTrackerMuonCands').at(0).hasPathName('HLT_Mu0_TkMu0_OST_Jpsi_Tight_v2') || "+
                                            "  triggerObjectMatchesByCollection('hltMuTkMuJpsiTrackerMuonCands').at(0).hasPathName('HLT_Mu0_TkMu0_OST_Jpsi_Tight_v3') )"),
)

LowPtTriggerFlagsEfficienciesTag = cms.PSet(

   ## Mu + Track

   Mu0_Track0_Jpsi_MU = cms.string("!triggerObjectMatchesByCollection('hltL3MuonCandidates').empty() && "+
                                   " triggerObjectMatchesByCollection('hltL3MuonCandidates').at(0).hasFilterLabel('hltMu0TrackJpsiTrackMassFiltered')"),
   Mu3_Track0_Jpsi_MU = cms.string("!triggerObjectMatchesByCollection('hltL3MuonCandidates').empty() && "+
                                   " triggerObjectMatchesByCollection('hltL3MuonCandidates').at(0).hasFilterLabel('hltMu3TrackJpsiTrackMassFiltered')"),
   Mu5_Track0_Jpsi_MU = cms.string("!triggerObjectMatchesByCollection('hltL3MuonCandidates').empty() && "+
                                   " triggerObjectMatchesByCollection('hltL3MuonCandidates').at(0).hasFilterLabel('hltMu5TrackJpsiTrackMassFiltered')"),
   Mu3_Track3_Jpsi_MU = cms.string("!triggerObjectMatchesByCollection('hltL3MuonCandidates').empty() && "+
                                   " triggerObjectMatchesByCollection('hltL3MuonCandidates').at(0).hasFilterLabel('hltMu3Track3JpsiTrackMassFiltered')"),
   Mu3_Track5_Jpsi_MU = cms.string("!triggerObjectMatchesByCollection('hltL3MuonCandidates').empty() && "+
                                   " triggerObjectMatchesByCollection('hltL3MuonCandidates').at(0).hasFilterLabel('hltMu3Track5JpsiTrackMassFiltered')"),
   ## Mu + L2Mu

   Mu0_L2Mu0_MU = cms.string("!triggerObjectMatchesByCollection('hltL3MuonCandidates').empty() && "+
                             " triggerObjectMatchesByCollection('hltL3MuonCandidates').at(0).hasFilterLabel('hltMu0L2Mu0L3Filtered0')"),
   Mu3_L2Mu0_MU = cms.string("!triggerObjectMatchesByCollection('hltL3MuonCandidates').empty() && "+
                             " triggerObjectMatchesByCollection('hltL3MuonCandidates').at(0).hasFilterLabel('hltMu3L2Mu0L3Filtered3')"),
   Mu5_L2Mu0_MU = cms.string("!triggerObjectMatchesByCollection('hltL3MuonCandidates').empty() && "+
                             " triggerObjectMatchesByCollection('hltL3MuonCandidates').at(0).hasFilterLabel('hltMu5L2Mu0L3Filtered5')"),

)


LowPtTriggerFlagsEfficienciesProbe = cms.PSet(


   ## Mu + Track

   Mu0_Track0_Jpsi_TK = cms.string("!triggerObjectMatchesByCollection('hltMuTrackJpsiCtfTrackCands').empty() && "+
                                   " triggerObjectMatchesByCollection('hltMuTrackJpsiCtfTrackCands').at(0).hasFilterLabel('hltMu0TrackJpsiTrackMassFiltered')"),
   Mu3_Track0_Jpsi_TK = cms.string("!triggerObjectMatchesByCollection('hltMuTrackJpsiCtfTrackCands').empty() && "+
                                   " triggerObjectMatchesByCollection('hltMuTrackJpsiCtfTrackCands').at(0).hasFilterLabel('hltMu3TrackJpsiTrackMassFiltered')"),
   Mu5_Track0_Jpsi_TK = cms.string("!triggerObjectMatchesByCollection('hltMuTrackJpsiCtfTrackCands').empty() && "+
                                   " triggerObjectMatchesByCollection('hltMuTrackJpsiCtfTrackCands').at(0).hasFilterLabel('hltMu5TrackJpsiTrackMassFiltered')"),
   Mu3_Track5_Jpsi_TK = cms.string("!triggerObjectMatchesByCollection('hltMuTrackJpsiCtfTrackCands').empty() && "+
                                   " triggerObjectMatchesByCollection('hltMuTrackJpsiCtfTrackCands').at(0).hasFilterLabel('hltMu3Track5JpsiTrackMassFiltered')"),
   Mu3_Track3_Jpsi_TK = cms.string("!triggerObjectMatchesByCollection('hltMuTrackJpsiCtfTrackCands').empty() && "+
                                   " triggerObjectMatchesByCollection('hltMuTrackJpsiCtfTrackCands').at(0).hasFilterLabel('hltMu3Track3JpsiTrackMassFiltered')"),

   MuX_Track0_Jpsi_Or_TK = cms.string("!triggerObjectMatchesByCollection('hltMuTrackJpsiCtfTrackCands').empty() && "+
                                      "(triggerObjectMatchesByCollection('hltMuTrackJpsiCtfTrackCands').at(0).hasFilterLabel('hltMu0TrackJpsiTrackMassFiltered')"+
                                      "||triggerObjectMatchesByCollection('hltMuTrackJpsiCtfTrackCands').at(0).hasFilterLabel('hltMu3TrackJpsiTrackMassFiltered')"+                                        
                                      "||triggerObjectMatchesByCollection('hltMuTrackJpsiCtfTrackCands').at(0).hasFilterLabel('hltMu5TrackJpsiTrackMassFiltered'))"
                                      ), 
   
   ## Mu + L2Mu
   MuX_L2Mu0_L2 = cms.string("!triggerObjectMatchesByCollection('hltL2MuonCandidates').empty() && "+
                             " triggerObjectMatchesByCollection('hltL2MuonCandidates').at(0).hasFilterLabel('hltDiMuonL2PreFiltered0')"),
)

LowPtTriggerFlagsEfficiencies = cms.PSet(LowPtTriggerFlagsEfficienciesTag,LowPtTriggerFlagsEfficienciesProbe)

def mkUnbiasCut(pset, triggers):
    clist = []
    for t in triggers:
        if hasattr(pset,t+"_MU"): clist.append(getattr(pset,t+"_MU").value())
        else: clist.append(getattr(pset,t).value())
    cut = " || " . join(["(%s)" % c for c in clist])
    return cut

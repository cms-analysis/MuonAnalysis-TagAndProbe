import FWCore.ParameterSet.Config as cms

process = cms.Process("Skim")

##     ____             __ _                       _     _           
##    / ___|___  _ __  / _(_) __ _ _   _ _ __ __ _| |__ | | ___  ___ 
##   | |   / _ \| '_ \| |_| |/ _` | | | | '__/ _` | '_ \| |/ _ \/ __|
##   | |__| (_) | | | |  _| | (_| | |_| | | | (_| | |_) | |  __/\__ \
##    \____\___/|_| |_|_| |_|\__, |\__,_|_|  \__,_|_.__/|_|\___||___/
##                           |___/                                   
##   
ptMin  = 1.5; # minimum pT to consider tracks and muons
ptTag  = 3.0; # minimum pT for the tag
etaTag = 2.5; # maximum |eta| for the tag
triggerProcess = 'HLT8E29'  # when running on 3.1.X MC
#triggerProcess = 'HLT'     # when testing on RelVals
triggerPath    = 'HLT_Mu3'  # Trigger name
massRangeMu  = (8.0, 12.0)
massRangeTk  = (8.0, 12.0) 
massRangeSta = (5.0, 14.0)



process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.source = cms.Source("PoolSource",  
    fileNames = cms.untracked.vstring(
        '/store/mc/Summer09/Upsilon1S/AODSIM/MC_31X_V3_AODSIM-v1/0011/D0A65840-A989-DE11-BC70-003048593AA6.root'
    )
)
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
## === GlobalTag is needed to read L1 infomration ====
#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#process.GlobalTag.globaltag = "MC_31X_V3::All"

##     ____  _ _                     _                  _    ___  ____  
##    / ___|| (_)_ __ ___  _ __ ___ (_)_ __   __ _     / \  / _ \|  _ \ 
##    \___ \| | | '_ ` _ \| '_ ` _ \| | '_ \ / _` |   / _ \| | | | | | |
##     ___) | | | | | | | | | | | | | | | | | (_| |  / ___ \ |_| | |_| |
##    |____/|_|_|_| |_| |_|_| |_| |_|_|_| |_|\__, | /_/   \_\___/|____/ 
##                                           |___/                      
## ==== pruner of GenParticles ====
process.genMuons = cms.EDProducer("GenParticlePruner",
    src = cms.InputTag("genParticles"),
    select = cms.vstring(
        "drop  *  ",                     # this is the default
        "++keep abs(pdgId) = 13",        # keep muons and their parents
        "drop pdgId == 21 && status = 2" # remove intermediate qcd spam carrying no flavour info
    )
)
## ==== Just ignore all the too low pt stuff ====
process.goodTracks = cms.EDFilter("TrackSelector",
    src = cms.InputTag("generalTracks"),
    cut = cms.string("pt > %f" % (ptMin,) ),
)

process.slimAOD = cms.Sequence(
    process.genMuons +
    process.goodTracks 
)



##    __  __                        ____      _       __  __                       
##   |  \/  | ___ _ __ __ _  ___   / ___|__ _| | ___ |  \/  |_   _  ___  _ __  ___ 
##   | |\/| |/ _ \ '__/ _` |/ _ \ | |   / _` | |/ _ \| |\/| | | | |/ _ \| '_ \/ __|
##   | |  | |  __/ | | (_| |  __/ | |__| (_| | | (_) | |  | | |_| | (_) | | | \__ \
##   |_|  |_|\___|_|  \__, |\___|  \____\__,_|_|\___/|_|  |_|\__,_|\___/|_| |_|___/
##                    |___/                                                        
##   
## ==== Merge CaloMuons into the collection of reco::Muons ====
from RecoMuon.MuonIdentification.calomuons_cfi import calomuons;
process.muons = cms.EDProducer("CaloMuonMerger",
    muons     = cms.InputTag("muons"), # half-dirty thing. it works as long as we're the first module using muons in the path
    caloMuons = cms.InputTag("calomuons"),
    minCaloCompatibility = calomuons.minCaloCompatibility
)

##    __  __       _          ____   _  _____   __  __                       
##   |  \/  | __ _| | _____  |  _ \ / \|_   _| |  \/  |_   _  ___  _ __  ___ 
##   | |\/| |/ _` | |/ / _ \ | |_) / _ \ | |   | |\/| | | | |/ _ \| '_ \/ __|
##   | |  | | (_| |   <  __/ |  __/ ___ \| |   | |  | | |_| | (_) | | | \__ \
##   |_|  |_|\__,_|_|\_\___| |_| /_/   \_\_|   |_|  |_|\__,_|\___/|_| |_|___/
##                                                                           
##   
### ==== Make PAT Muons ====
from PhysicsTools.PatAlgos.producersLayer1.muonProducer_cfi import allLayer1Muons
process.patMuonsWithoutTrigger = allLayer1Muons.clone(
    # embed the tracks, so we don't have to carry them around
    embedTrack          = True,
    embedCombinedMuon   = True,
    embedStandAloneMuon = True,
    # then switch off some features we don't need
    #addTeVRefits = False, ## <<--- this doesn't work. PAT bug ??
    embedPickyMuon = False,
    embedTpfmsMuon = False, 
    isolation = cms.PSet(),   # no extra isolation beyond what's in reco::Muon itself
    isoDeposits = cms.PSet(), # no heavy isodeposits
    addGenMatch = False,      # no mc: T&P doesn't take it from here anyway.
)
### ==== Unpack trigger, and match ====
process.load("PhysicsTools.PatAlgos.triggerLayer1.triggerProducer_cfi")
process.patTrigger.onlyStandAlone = True
process.patTrigger.processName    = triggerProcess
from PhysicsTools.PatAlgos.triggerLayer1.triggerMatcher_cfi import muonTriggerMatchHLT1MuonIso
process.muonMatchHLTMuX = muonTriggerMatchHLT1MuonIso.clone(
    src = 'patMuonsWithoutTrigger',
    pathNames = [ triggerPath ]
)
### ==== Embed ====
process.patMuonsWithTrigger = cms.EDProducer( "PATTriggerMatchMuonEmbedder",
    src     = cms.InputTag( "patMuonsWithoutTrigger" ),
    matches = cms.VInputTag( "muonMatchHLTMuX" ),
)
### ==== Select ====
process.patMuons = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("patMuonsWithTrigger"),
    cut = cms.string("pt > %f" % (ptMin,)), 
)
### ==== Sequence ====
process.patMuonSequence = cms.Sequence( 
    process.patMuonsWithoutTrigger *
    process.patTrigger * process.muonMatchHLTMuX * process.patMuonsWithTrigger *
    process.patMuons  
)

##    __  __       _                _   _                                                 
##   |  \/  | __ _| | _____    ___ | |_| |__   ___ _ __   _ __ ___  _   _  ___  _ __  ___ 
##   | |\/| |/ _` | |/ / _ \  / _ \| __| '_ \ / _ \ '__| | '_ ` _ \| | | |/ _ \| '_ \/ __|
##   | |  | | (_| |   <  __/ | (_) | |_| | | |  __/ |    | | | | | | |_| | (_) | | | \__ \
##   |_|  |_|\__,_|_|\_\___|  \___/ \__|_| |_|\___|_|    |_| |_| |_|\__,_|\___/|_| |_|___/
##                                                                                        
##   
## ==== For bare tracks, make candidates assuming the muon mass hypothesis ====
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi");
process.tkTracks  = cms.EDProducer("ConcreteChargedCandidateProducer", 
    src  = cms.InputTag("goodTracks"),      
    particleType = cms.string("mu+"),
) 
## ==== For skimming with standalone muons, use the raw standalone track, so that it doesn't take the tracker P4 if it's also a global or tracker muon ====
process.staTracks = cms.EDProducer("ConcreteChargedCandidateProducer", 
    src  = cms.InputTag("standAloneMuons","UpdatedAtVtx"), 
    particleType = cms.string("mu+"),
)
## ==== Muons that are not only standalone (for skimming only) ====
process.nonStaOnlyMuon = cms.EDFilter("PATMuonRefSelector",
    src = cms.InputTag("patMuons"),
    cut = cms.string("track.isNonnull"),
)

## ==== Golden muons to be used for tags. use PAT ones, so I can check HLT =====
PASS_HLT = "!triggerObjectMatchesByPath('%s').empty()" % (triggerPath,);
process.goldenMuons = cms.EDFilter("PATMuonRefSelector",
    src = cms.InputTag("patMuons"),
    cut = cms.string("isGlobalMuon && pt > %f && abs(eta) < %f && %s" % (ptTag, etaTag, PASS_HLT) ), 
)

process.otherMuons = cms.Sequence(
    process.nonStaOnlyMuon +
    process.goldenMuons +
    process.tkTracks +
    process.staTracks
)


process.allMuons = cms.Sequence(
    process.muons *  # this merges the CaloMuons
    process.patMuonSequence *
    process.otherMuons
)

##    __  __       _              _   ______      _ _     
##   |  \/  | __ _| | _____      | | / /  _ \ ___(_| )___ 
##   | |\/| |/ _` | |/ / _ \  _  | |/ /| |_) / __| |// __|
##   | |  | | (_| |   <  __/ | |_| / / |  __/\__ \ | \__ \
##   |_|  |_|\__,_|_|\_\___|  \___/_/  |_|   |___/_| |___/
##                                                        
##   
process.upsilonMu  = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("goldenMuons@+ nonStaOnlyMuon@-"),
    cut = cms.string("%f < mass < %f" % massRangeMu),
)
process.upsilonTk  = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("goldenMuons@+ tkTracks@-"),
    cut = cms.string("%f < mass < %f" % massRangeTk),
)
process.upsilonSta = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("goldenMuons@+ staTracks@-"),
    cut = cms.string("%f < mass < %f" % massRangeSta),
)

##    __  __       _                     _   _     
##   |  \/  | __ _(_)_ __    _ __   __ _| |_| |__  
##   | |\/| |/ _` | | '_ \  | '_ \ / _` | __| '_ \ 
##   | |  | | (_| | | | | | | |_) | (_| | |_| | | |
##   |_|  |_|\__,_|_|_| |_| | .__/ \__,_|\__|_| |_|
##                          |_|                    
##   
process.main = cms.Path(
    process.slimAOD  *    
    process.allMuons
)

##    ____  _    _             ____       _   _         
##   / ___|| | _(_)_ __ ___   |  _ \ __ _| |_| |__  ___ 
##   \___ \| |/ / | '_ ` _ \  | |_) / _` | __| '_ \/ __|
##    ___) |   <| | | | | | | |  __/ (_| | |_| | | \__ \
##   |____/|_|\_\_|_| |_| |_| |_|   \__,_|\__|_| |_|___/
##                                                      
##   
process.upsilonMuFilter  = cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag("upsilonMu"),
    minNumber = cms.uint32(1),
)
process.upsilonTkFilter  = process.upsilonMuFilter.clone(src = 'upsilonTk')
process.upsilonStaFilter = process.upsilonMuFilter.clone(src = 'upsilonSta')

process.Skim_upsilonMu  = cms.Path(process.upsilonMu  * process.upsilonMuFilter )
process.Skim_upsilonTk  = cms.Path(process.upsilonTk  * process.upsilonTkFilter )
process.Skim_upsilonSta = cms.Path(process.upsilonSta * process.upsilonStaFilter)

# I define a few other paths that are not enabled in the OutputModule but which are useful to see what is the
# efficiency of the various steps
# to run them and see the efficiency, uncomment them in the cms.Schedule
from HLTrigger.HLTfilters.hltHighLevel_cfi import hltHighLevel
process.filterHLTMuX = hltHighLevel.clone(
    HLTPaths = [ triggerPath  ],
    TriggerResultsTag = cms.InputTag("TriggerResults","",triggerProcess),
)
process.oneTagFilter  = process.upsilonMuFilter.clone(src = 'goldenMuons')
process.acceptanceFilter = cms.EDFilter("CandViewRefSelector",
    src = cms.InputTag("genMuons"),
    cut = cms.string(("pdgId == 553 && abs(daughter(0).pdgId) == 13 && ")+
                     ("abs(daughter(0).eta) < %f && abs(daughter(1).eta) < %f && " % (etaTag, etaTag))+
                     ("daughter(0).pt > %f && daughter(1).pt > %f && max(daughter(0).pt, daughter(1).pt) > %f" % (ptMin, ptMin, ptTag)) ),
    filter = cms.bool(True),
)
process.Check_OneTag  = cms.Path(process.filterHLTMuX * process.oneTagFilter)
# i clone these, as now they'll have different filters before...
process.filterHLTAcc  = process.filterHLTMuX.clone()
process.oneTagFilterAcc = process.oneTagFilter.clone()
process.upsilonMuFilterAcc = process.upsilonMuFilter.clone()
process.upsilonTkFilterAcc = process.upsilonTkFilter.clone()
process.upsilonStaFilterAcc = process.upsilonStaFilter.clone()
process.Check_Acceptance_Mu  = cms.Path(process.acceptanceFilter * process.filterHLTAcc * process.oneTagFilterAcc * process.upsilonMuFilterAcc)
process.Check_Acceptance_Tk  = cms.Path(process.acceptanceFilter * process.filterHLTAcc * process.oneTagFilterAcc * process.upsilonTkFilterAcc)
process.Check_Acceptance_Sta = cms.Path(process.acceptanceFilter * process.filterHLTAcc * process.oneTagFilterAcc * process.upsilonStaFilterAcc)


##     ___        _               _   
##    / _ \ _   _| |_ _ __  _   _| |_ 
##   | | | | | | | __| '_ \| | | | __|
##   | |_| | |_| | |_| |_) | |_| | |_ 
##    \___/ \__,_|\__| .__/ \__,_|\__|
##                   |_|              
##   
process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string("skimUpsilon.root"),
    outputCommands = cms.untracked.vstring(
        "drop *",
        "keep *_patMuons_*_Skim",
        "keep *_goodTracks_*_Skim",
        "keep *_genMuons_*_Skim",
        "keep recoTrackExtras_standAloneMuons_*_*",       ## track states at the muon system, used both by patMuons and standAloneMuons
        "keep recoTracks_standAloneMuons_UpdatedAtVtx_*", ## bare standalone muon tracks, using standalone muon momentum (with BS constraint)
        "keep edmTriggerResults_*_*_Skim",          ## to know which kind of skim channel got us the event   
       #"keep l1extraL1MuonParticles_l1extraParticles_*_*",  ## if we ever want to do L1 efficiency too ## <<--- Not in 3.1.X AODSIM
       # "keep *_upsilonMu_*_Skim", "keep *_upsilonTk_*_Skim", "keep *_upsilonSta_*_Skim",                       ## <<--- keep these for monitoring
    ),
    SelectEvents = cms.untracked.PSet( SelectEvents = cms.vstring(
        "Skim_upsilonMu",
        "Skim_upsilonTk",
        "Skim_upsilonSta",
    )),
)
process.end = cms.EndPath(process.out)

##    _____ _       _     _     _               _                   _               
##   |  ___(_)_ __ (_)___| |__ (_)_ __   __ _  | |_ ___  _   _  ___| |__   ___  ___ 
##   | |_  | | '_ \| / __| '_ \| | '_ \ / _` | | __/ _ \| | | |/ __| '_ \ / _ \/ __|
##   |  _| | | | | | \__ \ | | | | | | | (_| | | || (_) | |_| | (__| | | |  __/\__ \
##   |_|   |_|_| |_|_|___/_| |_|_|_| |_|\__, |  \__\___/ \__,_|\___|_| |_|\___||___/
##                                      |___/                                       
##   
process.schedule = cms.Schedule(
    process.main,
    process.Skim_upsilonMu,
    process.Skim_upsilonTk,
    process.Skim_upsilonSta,
    process.Check_OneTag, 
    process.Check_Acceptance_Mu, 
    process.Check_Acceptance_Tk, 
    process.Check_Acceptance_Sta, 
    process.end
)
#process.out.fileName = "/tmp/gpetrucc/skimUpsilon.root"
process.maxEvents.input = 10000

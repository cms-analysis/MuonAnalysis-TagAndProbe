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
massRangeMu  = (2.0, 4.0)
massRangeTk  = (2.0, 4.0) 
massRangeSta = (2.0, 5.0)



process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.source = cms.Source("PoolSource",  
    fileNames = cms.untracked.vstring(
        'file:/afs/cern.ch/user/g/gpetrucc/scratch0/huntForRedOctober/CMSSW_3_1_2/src/MuonAnalysis/TagAndProbe/test/crab/ppMuX_AODSIM.root'
    )
)
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
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
process.nonStaOnlyMuon = cms.EDFilter("MuonRefSelector",
    src = cms.InputTag("muons"),
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
process.jpsiMu  = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("goldenMuons@+ nonStaOnlyMuon@-"),
    cut = cms.string("%f < mass < %f" % massRangeMu),
)
process.jpsiTk  = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("goldenMuons@+ tkTracks@-"),
    cut = cms.string("%f < mass < %f" % massRangeTk),
)
process.jpsiSta = cms.EDProducer("CandViewShallowCloneCombiner",
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
process.jpsiMuFilter  = cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag("jpsiMu"),
    minNumber = cms.uint32(1),
)
process.jpsiTkFilter  = process.jpsiMuFilter.clone(src = 'jpsiTk')
process.jpsiStaFilter = process.jpsiMuFilter.clone(src = 'jpsiSta')

process.Skim_jpsiMu  = cms.Path(process.jpsiMu  * process.jpsiMuFilter )
process.Skim_jpsiTk  = cms.Path(process.jpsiTk  * process.jpsiTkFilter )
process.Skim_jpsiSta = cms.Path(process.jpsiSta * process.jpsiStaFilter)

# I define also a path that just skims on the trigger bit, to see what's the efficiency.
# it won't be enabled in the PoolOutputModule
from HLTrigger.HLTfilters.hltHighLevel_cfi import hltHighLevel
process.filterHLTMuX = hltHighLevel.clone(
    HLTPaths = [ triggerPath  ],
    TriggerResultsTag = cms.InputTag("TriggerResults","",triggerProcess),
)
process.Skim_HLT_MuX  = cms.Path(process.filterHLTMuX)


##     ___        _               _   
##    / _ \ _   _| |_ _ __  _   _| |_ 
##   | | | | | | | __| '_ \| | | | __|
##   | |_| | |_| | |_| |_) | |_| | |_ 
##    \___/ \__,_|\__| .__/ \__,_|\__|
##                   |_|              
##   
process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string("skimJPsi.root"),
    outputCommands = cms.untracked.vstring(
        "drop *",
        "keep *_patMuons_*_Skim",
        "keep *_goodTracks_*_Skim",
        "keep *_genMuons_*_Skim",
        "keep recoTrackExtras_standAloneMuons_*_*", ## track states at the muon system
        "keep edmTriggerResults_*_*_Skim",          ## to know which kind of skim channel got us the event   
        "keep l1extraL1MuonParticles_hltL1extraParticles_*_*", ## if we ever want to do L1 efficiency too
    ),
    SelectEvents = cms.untracked.PSet( SelectEvents = cms.vstring(
        "Skim_jpsiMu",
        "Skim_jpsiTk",
        "Skim_jpsiSta",
    )),
    dropMetaData = cms.untracked.string("PRIOR"), ## We don't need to keep the metadata for the previous processes, 
                                                  ## as we can get the cfgs from DBS
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
    process.Skim_jpsiMu,
    process.Skim_jpsiTk,
    process.Skim_jpsiSta,
    process.Skim_HLT_MuX, 
    process.end
)
#process.out.fileName = "/tmp/gpetrucc/skimJPsi.root"
#process.maxEvents.input = 5000

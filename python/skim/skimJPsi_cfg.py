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
#triggerProcess = 'HLT8E29'  # when running on 3.1.X MC
triggerProcess = 'HLT'       # when  running on 3.3.X MC or RelVals
triggerPath    = 'HLT_Mu3'   # Trigger name
massRangeMu  = (2.0, 4.0)
massRangeTk  = (2.0, 4.0) 
massRangeSta = (2.0, 5.0)



process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.source = cms.Source("PoolSource",  
    fileNames = cms.untracked.vstring(
        '/store/relval/CMSSW_3_4_1/RelValJpsiMM/GEN-SIM-RECO/STARTUP3X_V14-v1/0004/DC4D9E50-82ED-DE11-BF03-003048D375AA.root',
        '/store/relval/CMSSW_3_4_1/RelValJpsiMM/GEN-SIM-RECO/STARTUP3X_V14-v1/0004/88F9B6DC-85ED-DE11-A295-001D09F2514F.root',
        '/store/relval/CMSSW_3_4_1/RelValJpsiMM/GEN-SIM-RECO/STARTUP3X_V14-v1/0004/66B9D232-85ED-DE11-B5EC-0019B9F7312C.root',
        '/store/relval/CMSSW_3_4_1/RelValJpsiMM/GEN-SIM-RECO/STARTUP3X_V14-v1/0004/509FE9BB-83ED-DE11-A25F-001D09F253C0.root',
        '/store/relval/CMSSW_3_4_1/RelValJpsiMM/GEN-SIM-RECO/STARTUP3X_V14-v1/0004/4294CBAC-B5ED-DE11-B7E8-003048D2BE08.root',
    )
)
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
## === GlobalTag is needed to read L1 infomration ====
#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#process.GlobalTag.globaltag = "STARTUP3X_V15::All"

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


##    __  __       _          ____   _  _____   __  __                       
##   |  \/  | __ _| | _____  |  _ \ / \|_   _| |  \/  |_   _  ___  _ __  ___ 
##   | |\/| |/ _` | |/ / _ \ | |_) / _ \ | |   | |\/| | | | |/ _ \| '_ \/ __|
##   | |  | | (_| |   <  __/ |  __/ ___ \| |   | |  | | |_| | (_) | | | \__ \
##   |_|  |_|\__,_|_|\_\___| |_| /_/   \_\_|   |_|  |_|\__,_|\___/|_| |_|___/
##
##                                                                           
##  We take them from Onia2MuMuPAT but we make a few changes
##
process.load("HeavyFlavorAnalysis.Onia2MuMu.onia2MuMuPAT_cff");
## put a PT cut on the muons
process.patMuons.cut = "pt > %f" % (ptMin,);
## use the genMuons as MC source, so that we keep them and have the correct mother refs
process.muonMatch.matched = 'genMuons'
## switch off genParticle embedding, as we keep the genMuons collection
## also switch off standalone muon track embedding, as we keep it separately
process.patMuonsWithoutTrigger.embedGenMatch = False
process.patMuonsWithoutTrigger.embedStandAloneMuon = False

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
process.allJPsis = cms.Sequence(
    process.jpsiMu + process.jpsiTk + process.jpsiSta
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
    process.allMuons *
    process.allJPsis
)

##     ____      _ _ _     _               _____                 _     ____       _           _   _             
##    / ___|___ | | (_)___(_) ___  _ __   | ____|_   _____ _ __ | |_  / ___|  ___| | ___  ___| |_(_) ___  _ __  
##   | |   / _ \| | | / __| |/ _ \| '_ \  |  _| \ \ / / _ \ '_ \| __| \___ \ / _ \ |/ _ \/ __| __| |/ _ \| '_ \ 
##   | |__| (_) | | | \__ \ | (_) | | | | | |___ \ V /  __/ | | | |_   ___) |  __/ |  __/ (__| |_| | (_) | | | |
##    \____\___/|_|_|_|___/_|\___/|_| |_| |_____| \_/ \___|_| |_|\__| |____/ \___|_|\___|\___|\__|_|\___/|_| |_|
##                                                                                                              
##   

# Some event filters
process.load('L1TriggerConfig.L1GtConfigProducers.L1GtTriggerMaskTechTrigConfig_cff')
from HLTrigger.HLTfilters.hltLevel1GTSeed_cfi import hltLevel1GTSeed
process.bit40    = hltLevel1GTSeed.clone(L1TechTriggerSeeding = cms.bool(True), L1SeedsLogicalExpression = cms.string('(40 OR 41)'))
process.haloVeto = hltLevel1GTSeed.clone(L1TechTriggerSeeding = cms.bool(True), L1SeedsLogicalExpression = cms.string('NOT (36 OR 37 OR 38 OR 39)'))
process.bptxAnd  = hltLevel1GTSeed.clone(L1TechTriggerSeeding = cms.bool(True), L1SeedsLogicalExpression = cms.string('0'))

from HLTrigger.HLTfilters.hltHighLevelDev_cfi import hltHighLevelDev
process.physDecl = hltHighLevelDev.clone(HLTPaths = ['HLT_PhysicsDeclared'], HLTPathsPrescales = [1])

process.oneGoodVertexFilter = cms.EDFilter("VertexSelector",
   src = cms.InputTag("offlinePrimaryVertices"),
   cut = cms.string("!isFake && tracksSize > 3 && abs(z) <= 15 && position.Rho <= 2"),
   filter = cms.bool(True),   # otherwise it won't filter the events, just produce an empty vertex collection.
)

process.noScraping = cms.EDFilter("FilterOutScraping",
    applyfilter = cms.untracked.bool(True),
    debugOn = cms.untracked.bool(False), ## Or 'True' to get some per-event info
    numtrack = cms.untracked.uint32(10),
    thresh = cms.untracked.double(0.2)
)

#### Current configuration
## rationale:
##   - bptxAnd and haloVeto are needed only on data (MC has no halo nor beam-gas)
##   - oneGoodVertexFilter should also cut away cosmics in coincidence with BPTX
##   - noScraping shouldn't cut anything, but we keep it for safety
process.collisionFilterMC   = cms.Sequence(process.oneGoodVertexFilter * process.noScraping)
process.collisionFilterData = cms.Sequence(process.bptxAnd * process.haloVeto * process.oneGoodVertexFilter * process.noScraping)
process.collisionFilter     = cms.Sequence(process.collisionFilterMC)

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

process.Skim_jpsiMu  = cms.Path(process.collisionFilter + process.jpsiMuFilter )
process.Skim_jpsiTk  = cms.Path(process.collisionFilter + process.jpsiTkFilter )
process.Skim_jpsiSta = cms.Path(process.collisionFilter + process.jpsiStaFilter)

# I define a few other paths that are not enabled in the OutputModule but which are useful to see what is the
# efficiency of the various steps
# to run them and see the efficiency, uncomment them in the cms.Schedule
from HLTrigger.HLTfilters.hltHighLevel_cfi import hltHighLevel
process.filterHLTMuX = hltHighLevel.clone(
    HLTPaths = [ triggerPath  ],
    TriggerResultsTag = cms.InputTag("TriggerResults","",triggerProcess),
)
process.oneTagFilter  = process.jpsiMuFilter.clone(src = 'goldenMuons')

process.Check_Collisions = cms.Path(process.collisionFilter)
process.Check_OneTag     = cms.Path(process.collisionFilter + process.filterHLTMuX + process.oneTagFilter)

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
        "keep recoTrackExtras_standAloneMuons_*_*",          ## track states at the muon system, used both by patMuons and standAloneMuons
        "keep recoTracks_standAloneMuons_UpdatedAtVtx_*",    ## bare standalone muon tracks, using standalone muon momentum (with BS constraint)
        "keep edmTriggerResults_*_*_Skim",                   ## to know which kind of skim channel got us the event   
        "keep l1extraL1MuonParticles_l1extraParticles_*_*",  ## if we ever want to do L1 efficiency too ## <<--- Not in 3.1.X AODSIM
        "keep *_offlinePrimaryVertices__*",                  ## vertices and BS are not very useful on MC
        "keep *_offlineBeamSpot__*",                         ## but they can be important on data
       # "keep *_jpsiMu_*_Skim", "keep *_jpsiTk_*_Skim", "keep *_jpsiSta_*_Skim",                       ## <<--- keep these for monitoring
    ),
    SelectEvents = cms.untracked.PSet( SelectEvents = cms.vstring(
        "Skim_jpsiMu",
        "Skim_jpsiTk",
        "Skim_jpsiSta",
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
    process.Skim_jpsiMu,
    process.Skim_jpsiTk,
    process.Skim_jpsiSta,
    process.Check_Collisions, 
    process.Check_OneTag, 
    process.end
)



import FWCore.ParameterSet.Config as cms

ptMinCut       = 'pt > 2 || (abs(eta) > 1 && p > 2)';
triggerPath1Mu = 'HLT_Mu3'                 # Single-muon Trigger name
triggerFilt1Mu = 'hltSingleMu3L3Filtered3' # Single-muon Trigger last filter name
triggerPath2Mu = 'HLT_L1DoubleMuOpen'                  # Double-muon Trigger name
triggerFilt2Mu = 'hltDoubleMuLevel1PathL1OpenFiltered' # Double-muon Trigger name
massRangeMu  = (2.6, 3.6)
massRangeTk  = (2.6, 3.6) 
massRangeSta = (2.0, 5.0)


##     ____  _ _                     _                  _    ___  ____  
##    / ___|| (_)_ __ ___  _ __ ___ (_)_ __   __ _     / \  / _ \|  _ \ 
##    \___ \| | | '_ ` _ \| '_ ` _ \| | '_ \ / _` |   / _ \| | | | | | |
##     ___) | | | | | | | | | | | | | | | | | (_| |  / ___ \ |_| | |_| |
##    |____/|_|_|_| |_| |_|_| |_| |_|_|_| |_|\__, | /_/   \_\___/|____/ 
##                                           |___/                      
## ==== pruner of GenParticles ====
genMuons = cms.EDProducer("GenParticlePruner",
    src = cms.InputTag("genParticles"),
    select = cms.vstring(
        "drop  *  ",                     # this is the default
        "++keep abs(pdgId) = 13",        # keep muons and their parents
        "drop pdgId == 21 && status = 2" # remove intermediate qcd spam carrying no flavour info
    )
)
## ==== Just ignore all the too low pt stuff ====
goodTracks = cms.EDFilter("TrackSelector",
    src = cms.InputTag("generalTracks"),
    cut = cms.string(ptMinCut),
)

## ==== Slim PAT jets, for B analysis that need dr(mu,jet) ====
from PhysicsTools.PatAlgos.producersLayer1.jetProducer_cff import patJetCorrFactors
from PhysicsTools.PatAlgos.producersLayer1.jetProducer_cff import patJets as ak5CaloJetsPAT
patJetCorrFactors.corrSample = "Summer09_7TeV_ReReco332" ## << TO BE IMPROVED
patJetCorrFactors.corrLevels.L5Flavor = 'none' # save space
patJetCorrFactors.corrLevels.L7Parton = 'none' # save space
ak5CaloJetsPAT.addBTagInfo         = False # save
ak5CaloJetsPAT.addAssociatedTracks = False # more
ak5CaloJetsPAT.addJetCharge        = False # space
ak5CaloJetsPAT.addGenPartonMatch = False # save space
ak5CaloJetsPAT.addGenJetMatch    = False # and not
ak5CaloJetsPAT.getJetMCFlavour   = False # in data
patJets = cms.EDFilter("PATJetSelector",
    src = cms.InputTag("ak5CaloJetsPAT"),
    cut = cms.string("pt > 20 && abs(eta) < 3 && (emEnergyFraction() > 0.01 && jetID.n90Hits > 1 && jetID.fHPD < 0.98)"),
)

slimAOD = cms.Sequence(
    genMuons +
    goodTracks +
    (patJetCorrFactors * ak5CaloJetsPAT) * patJets
)


##    __  __       _          ____   _  _____   __  __                       
##   |  \/  | __ _| | _____  |  _ \ / \|_   _| |  \/  |_   _  ___  _ __  ___ 
##   | |\/| |/ _` | |/ / _ \ | |_) / _ \ | |   | |\/| | | | |/ _ \| '_ \/ __|
##   | |  | | (_| |   <  __/ |  __/ ___ \| |   | |  | | |_| | (_) | | | \__ \
##   |_|  |_|\__,_|_|\_\___| |_| /_/   \_\_|   |_|  |_|\__,_|\___/|_| |_|___/
##
##                                                                           
## ==== Merge CaloMuons into the collection of reco::Muons  ====
from RecoMuon.MuonIdentification.calomuons_cfi import calomuons;
mergedMuons = cms.EDProducer("CaloMuonMerger",
    muons     = cms.InputTag("muons"), 
    caloMuons = cms.InputTag("calomuons"),
    minCaloCompatibility = calomuons.minCaloCompatibility
)

## ==== PAT Muons with trigger matching ====
from MuonAnalysis.MuonAssociators.patMuonsWithTrigger_8E29_cff import *
## use merged collection
muonL1Info.src = 'mergedMuons'
patMuonsWithoutTrigger.muonSource = 'mergedMuons'
## also switch off standalone muon track embedding, as we keep it separately
patMuonsWithoutTrigger.embedStandAloneMuon = False
## also add matching with with di-muon triggers
patTriggerMatching.replace(patTriggerMatchers1Mu, patTriggerMatchers1Mu + patTriggerMatchers2Mu)
patMuonsWithTrigger.matches += patTriggerMatchers2MuInputTags
## finally make a collection of patMuons possibly applying some filter
patMuons = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("patMuonsWithTrigger"),
    cut = cms.string(ptMinCut),
)
patMuonSequence = cms.Sequence(
    mergedMuons *
    patMuonsWithTriggerSequence *
    patMuons
)

##    __  __       _                _   _                                                 
##   |  \/  | __ _| | _____    ___ | |_| |__   ___ _ __   _ __ ___  _   _  ___  _ __  ___ 
##   | |\/| |/ _` | |/ / _ \  / _ \| __| '_ \ / _ \ '__| | '_ ` _ \| | | |/ _ \| '_ \/ __|
##   | |  | | (_| |   <  __/ | (_) | |_| | | |  __/ |    | | | | | | |_| | (_) | | | \__ \
##   |_|  |_|\__,_|_|\_\___|  \___/ \__|_| |_|\___|_|    |_| |_| |_|\__,_|\___/|_| |_|___/
##                                                                                        
##   
## ==== For bare tracks, make candidates assuming the muon mass hypothesis ====
from SimGeneral.HepPDTESSource.pythiapdt_cfi import *
tkTracks  = cms.EDProducer("ConcreteChargedCandidateProducer", 
    src  = cms.InputTag("goodTracks"),      
    particleType = cms.string("mu+"),
) 
## ==== For skimming with standalone muons, use the raw standalone track, so that it doesn't take the tracker P4 if it's also a global or tracker muon ====
staTracks = cms.EDProducer("ConcreteChargedCandidateProducer", 
    src  = cms.InputTag("standAloneMuons","UpdatedAtVtx"), 
    particleType = cms.string("mu+"),
)
## ==== Muons that are not only standalone (for skimming only) ====
nonStaOnlyMuon = cms.EDFilter("PATMuonRefSelector",
    src = cms.InputTag("patMuons"),
    cut = cms.string("track.isNonnull"),
)

## ==== Golden muons to be used for tags. use PAT ones, so I can check HLT =====
PASS_HLT_1 = "!triggerObjectMatchesByFilter('%s').empty()" % (triggerFilt1Mu,);
PASS_HLT_2 = "!triggerObjectMatchesByFilter('%s').empty()" % (triggerFilt2Mu,);
tags1Mu = cms.EDFilter("PATMuonRefSelector",
    src = cms.InputTag("patMuons"),
    cut = cms.string("isGlobalMuon && "+ PASS_HLT_1), 
)
tags2Mu = cms.EDFilter("PATMuonRefSelector",
    src = cms.InputTag("patMuons"),
    cut = cms.string("isGlobalMuon && " + PASS_HLT_2), 
)

otherMuons = cms.Sequence(
    nonStaOnlyMuon +
    tags1Mu +
    tags2Mu +
    tkTracks +
    staTracks
)

allMuons = cms.Sequence(
    patMuonSequence *
    otherMuons
)

##    __  __       _              _   ______      _ _     
##   |  \/  | __ _| | _____      | | / /  _ \ ___(_| )___ 
##   | |\/| |/ _` | |/ / _ \  _  | |/ /| |_) / __| |// __|
##   | |  | | (_| |   <  __/ | |_| / / |  __/\__ \ | \__ \
##   |_|  |_|\__,_|_|\_\___|  \___/_/  |_|   |___/_| |___/
##                                                        
##   
jpsiMu  = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("tags1Mu@+ nonStaOnlyMuon@-"),
    cut = cms.string("%f < mass < %f" % massRangeMu),
)
jpsiTk  = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("tags1Mu@+ tkTracks@-"),
    cut = cms.string("%f < mass < %f" % massRangeTk),
)
jpsiSta = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("tags2Mu@+ staTracks@-"),
    cut = cms.string("%f < mass < %f" % massRangeSta),
)
allJPsis = cms.Sequence(
    jpsiMu + jpsiTk + jpsiSta
)

##    __  __       _                     _   _     
##   |  \/  | __ _(_)_ __    _ __   __ _| |_| |__  
##   | |\/| |/ _` | | '_ \  | '_ \ / _` | __| '_ \ 
##   | |  | | (_| | | | | | | |_) | (_| | |_| | | |
##   |_|  |_|\__,_|_|_| |_| | .__/ \__,_|\__|_| |_|
##                          |_|                    
##   
jpsiSkimMainSequence = cms.Sequence(
    slimAOD  *    
    allMuons *
    allJPsis
)

##    ____  _    _             ____       _   _         
##   / ___|| | _(_)_ __ ___   |  _ \ __ _| |_| |__  ___ 
##   \___ \| |/ / | '_ ` _ \  | |_) / _` | __| '_ \/ __|
##    ___) |   <| | | | | | | |  __/ (_| | |_| | | \__ \
##   |____/|_|\_\_|_| |_| |_| |_|   \__,_|\__|_| |_|___/
##                                                      
##   
jpsiMuFilter  = cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag("jpsiMu"),
    minNumber = cms.uint32(1),
)
jpsiTkFilter  = jpsiMuFilter.clone(src = 'jpsiTk')
jpsiStaFilter = jpsiMuFilter.clone(src = 'jpsiSta')


##     ___        _               _   
##    / _ \ _   _| |_ _ __  _   _| |_ 
##   | | | | | | | __| '_ \| | | | __|
##   | |_| | |_| | |_| |_) | |_| | |_ 
##    \___/ \__,_|\__| .__/ \__,_|\__|
##                   |_|              
##   
jpsiSkimOut = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string("skimJPsi.root"),
    #fileName = cms.untracked.string("/tmp/gpetrucc/skimJPsi_ppMuX.root"),
    outputCommands = cms.untracked.vstring(
        "drop *",
        "keep *_patMuons_*_Skim",
        "keep *_goodTracks_*_Skim",
        "keep *_genMuons_*_Skim",
        "keep recoTrackExtras_standAloneMuons_*_*",          ## track states at the muon system, used both by patMuons and standAloneMuons
        "keep recoTracks_standAloneMuons_UpdatedAtVtx_*",    ## bare standalone muon tracks, using standalone muon momentum (with BS constraint)
        "keep edmTriggerResults_*_*_Skim",                   ## to know which kind of skim channel got us the event   
        "keep edmTriggerResults_*_*_HLT",                    ## to know what got us on tape
        "keep l1extraL1MuonParticles_l1extraParticles_*_*",  ## if we ever want to do L1 efficiency too ## <<--- Not in 3.1.X AODSIM
        "keep *_offlinePrimaryVertices__*",                  ## vertices and BS are not very useful on MC
        "keep *_offlineBeamSpot__*",                         ## but they can be important on data
        "keep patJets_patJets__*",                           ## for inclusive B and for top backgrounds
       #"keep *_jpsiMu_*_Skim", "keep *_jpsiTk_*_Skim", "keep *_jpsiSta_*_Skim",                       ## <<--- keep these for monitoring
    ),
    SelectEvents = cms.untracked.PSet( SelectEvents = cms.vstring(
        "Skim_jpsiMu",
        "Skim_jpsiTk",
        "Skim_jpsiSta",
    )),
)

def changeTriggerProcessName(process, triggerProcessName, oldProcessName="HLT"):
    from MuonAnalysis.MuonAssociators.patMuonsWithTrigger_8E29_cff import changeTriggerProcessName as changeTriggerProcessNameBase
    changeTriggerProcessNameBase(process, triggerProcessName, oldProcessName)
    if hasattr(process, 'filterHLT1Mu'): process.filterHLT1Mu.TriggerResultsTag = cms.InputTag('TriggerResults','',triggerProcessName)
    if hasattr(process, 'filterHLT2Mu'): process.filterHLT2Mu.TriggerResultsTag = cms.InputTag('TriggerResults','',triggerProcessName)
    process.jpsiSkimOut.outputCommands += [
        "drop edmTriggerResults_*_*_"+oldProcessName,       ## to know what got us on tape
        "keep edmTriggerResults_*_*_"+triggerProcessName,   ## to know what got us on tape
    ]
    
def Summer09_Trigger(process):
    changeTriggerProcessName(process, 'HLT8E29', process.patTrigger.processName.value())

def Spring10ReDigi_Trigger(process):
    changeTriggerProcessName(process, 'HLT8E29', process.patTrigger.processName.value())


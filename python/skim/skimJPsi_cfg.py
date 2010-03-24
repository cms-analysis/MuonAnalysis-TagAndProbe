import FWCore.ParameterSet.Config as cms

process = cms.Process("Skim")


process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.source = cms.Source("PoolSource",  
    fileNames = cms.untracked.vstring(
	'/store/relval/CMSSW_3_5_4/RelValJpsiMM/GEN-SIM-RECO/START3X_V24-v1/0004/1244DD00-2D2C-DF11-8FA3-002618943977.root',
	'/store/relval/CMSSW_3_5_4/RelValJpsiMM/GEN-SIM-RECO/START3X_V24-v1/0004/1217EA49-A52B-DF11-9081-00304867BFAE.root',
	'/store/relval/CMSSW_3_5_4/RelValJpsiMM/GEN-SIM-RECO/START3X_V24-v1/0003/BE67D561-A02B-DF11-9FE7-001A928116B4.root',
	'/store/relval/CMSSW_3_5_4/RelValJpsiMM/GEN-SIM-RECO/START3X_V24-v1/0003/B892E028-9F2B-DF11-96DD-003048679164.root',
	'/store/relval/CMSSW_3_5_4/RelValJpsiMM/GEN-SIM-RECO/START3X_V24-v1/0003/5C692C2C-9E2B-DF11-B817-002618943869.root',
    )
)
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(5000))
## === conditions are needed to access geometry for L1 matching ====
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load("Configuration.StandardSequences.GeometryExtended_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.GlobalTag.globaltag = "MC_3XY_V25::All"  # 3.1.X MC is Ideal
#process.GlobalTag.globaltag = "START3X_V25::All"

process.load("MuonAnalysis.TagAndProbe.skim.skimJPsi_cff")
# we also do a regular python import to get some variables like the names of the triggers, which are not exported by process.load
from MuonAnalysis.TagAndProbe.skim.skimJPsi_cff import *

##    __  __       _                     _   _     
##   |  \/  | __ _(_)_ __    _ __   __ _| |_| |__  
##   | |\/| |/ _` | | '_ \  | '_ \ / _` | __| '_ \ 
##   | |  | | (_| | | | | | | |_) | (_| | |_| | | |
##   |_|  |_|\__,_|_|_| |_| | .__/ \__,_|\__|_| |_|
##                          |_|                    
##   
process.main = cms.Path( process.jpsiSkimMainSequence )

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
## Note: we don't use the collision filter when running on 'high luminosity' MC
#process.Skim_jpsiMu  = cms.Path(process.collisionFilter + process.jpsiMuFilter )
#process.Skim_jpsiTk  = cms.Path(process.collisionFilter + process.jpsiTkFilter )
#process.Skim_jpsiSta = cms.Path(process.collisionFilter + process.jpsiStaFilter)
process.Skim_jpsiMu  = cms.Path(process.jpsiMuFilter )
process.Skim_jpsiTk  = cms.Path(process.jpsiTkFilter )
process.Skim_jpsiSta = cms.Path(process.jpsiStaFilter)

# I define a few other paths that are not enabled in the OutputModule but which are useful to see what is the
# efficiency of the various steps
# to run them and see the efficiency, uncomment them in the cms.Schedule
from HLTrigger.HLTfilters.hltHighLevel_cfi import hltHighLevel
process.filterHLT1Mu = hltHighLevel.clone(
    HLTPaths = [ triggerPath1Mu  ],
    TriggerResultsTag = cms.InputTag("TriggerResults","","HLT"),
)

process.filterHLT2Mu = process.filterHLT1Mu.clone(HLTPaths = [ triggerPath2Mu ])
process.oneTag1MuFilter  = process.jpsiMuFilter.clone(src = 'tags1Mu')
process.oneTag2MuFilter  = process.jpsiMuFilter.clone(src = 'tags2Mu')


#process.Check_Collisions = cms.Path(process.collisionFilter)
process.Check_HLT1Mu  = cms.Path(process.filterHLT1Mu)
process.Check_HLT2Mu  = cms.Path(process.filterHLT2Mu)
process.Check_OneTag1Mu  = cms.Path(process.oneTag1MuFilter)
process.Check_OneTag2Mu  = cms.Path(process.oneTag2MuFilter)

##     ___        _               _   
##    / _ \ _   _| |_ _ __  _   _| |_ 
##   | | | | | | | __| '_ \| | | | __|
##   | |_| | |_| | |_| |_) | |_| | |_ 
##    \___/ \__,_|\__| .__/ \__,_|\__|
##                   |_|              
##   
process.end = cms.EndPath(process.jpsiSkimOut)

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
    #process.Check_Collisions, 
    process.Check_OneTag1Mu, 
    process.Check_OneTag2Mu, 
    process.Check_HLT1Mu, 
    process.Check_HLT2Mu, 
    process.end
)

##    _____         _       
##   |_   _|__  ___| |_ ___ 
##     | |/ _ \/ __| __/ __|
##     | |  __/\__ \ |_\__ \
##     |_|\___||___/\__|___/
##                          
##   

import sys
args = sys.argv[1:]
if (sys.argv[0] == "cmsRun"): args =sys.argv[2:]
test = "Summer09_Signal"
if len(args) > 0:
    test = args[0]
    print "Will run test '%s'" % (test,)

if test == "Summer09_Signal":
    Summer09_Trigger(process)
    process.source.fileNames = [ 'root://pcmssd12.cern.ch//data/mc/7TeV/Summer09/JPsiMuMu_AODSIM_72B0D83F-908B-DE11-9A88-001D09646131.root' ]
    process.jpsiSkimOut.fileName = 'JPsiMuMu_skimJPsi_Summmer09.root'
elif test == "Summer09_ppMuX":
    Summer09_Trigger(process)
    process.source.fileNames = [ 'root://pcmssd12.cern.ch//data/mc/7TeV/Summer09/ppMuX_AODSIM_EAB08D8C-78A6-DE11-88FD-001AA009562B.root' ]
    process.jpsiSkimOut.fileName = 'ppMuX_skimJPsi_Summmer09.root'
elif test == "RelVal354":
    process.jpsiSkimOut.fileName = 'JPsiMuMu_skimJPsi_RelVal354.root'
elif test == "OniaTrigger":
    process.source.fileNames = [ 'root://pcmssd12.cern.ch//data/gpetrucc/PYTHIA6_7TeV_JPsiWithFSR_MC_3XY_V16_1E31_HLT_12.root' ]
    process.source.duplicateCheckMode = cms.untracked.string("noDuplicateCheck")
    process.jpsiSkimOut.fileName = 'JPsiMuMu_skimJPsi_OniaTrigger.root'
else:
    raise ValueError, "Unknown test '%s'" % (test,)

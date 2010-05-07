import FWCore.ParameterSet.Config as cms

process = cms.Process("TagProbe")

process.load('FWCore.MessageService.MessageLogger_cfi')
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.source = cms.Source("PoolSource", 
    fileNames = cms.untracked.vstring(
        'rfio:/castor/cern.ch/user/g/gpetrucc/7TeV/DATA/DATA_skimJPsiLoose_fromApr20MuonSkim-v2.root',
        'rfio:/castor/cern.ch/user/g/gpetrucc/7TeV/DATA/DATA_skimJPsiLoose_fromMuonSkimV9_upToApr28-v2.root',
        #'root://pcmssd12.cern.ch//data/gpetrucc/7TeV/tnp/DATA_skimJPsiLoose_fromApr20MuonSkim-v2.root',
        #'root://pcmssd12.cern.ch//data/gpetrucc/7TeV/tnp/DATA_skimJPsiLoose_fromMuonSkimV9_upToApr28-v2.root',
    )
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )    

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Geometry.CommonDetUnit.globalTrackingGeometry_cfi")
process.load("TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorAny_cfi")
process.load("TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorAlong_cfi")
process.load("TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorOpposite_cfi")
process.GlobalTag.globaltag = cms.string('GR_R_35X_V8::All')

process.load("MuonAnalysis.TagAndProbe.jpsi.tp_from_skim_Tracking_cff")
process.load("MuonAnalysis.TagAndProbe.jpsi.tp_from_skim_MuonID_cff")
process.load("MuonAnalysis.TagAndProbe.jpsi.tp_from_skim_Trigger_cff")

process.tagAndProbe = cms.Path( 
    process.tnpCommonSequence    *
    ( process.tnpSequenceTracking  +
      process.tnpSequenceMuonID    +
      process.tnpSequenceMuonIDVtx +
      process.tnpSequenceTrigger   )
)

from MuonAnalysis.TagAndProbe.jpsi.tp_from_skim_common_cff import *

## Remove all MC matcher modules from the sequence
for K in process.filters_().keys() + process.producers_().keys():
    V = getattr(process,K) #; print K, V.type_()
    if V.type_() == "MCTruthDeltaRMatcherNew":
        process.tagAndProbe.remove(V)

## Set 'isMC' to False in all T&P tree producers
for K,V in allTPTreeProducers(process): V.isMC = False

## Redefine the tags requiring only L1SingleMuOpen
process.tagMuons1Mu.cut = PASSING_GLB_CUT + " && !triggerObjectMatchesByFilter('hltL1MuOpenL1Filtered0').empty() && " + TRACK_CUTS;
process.tagMuons2Mu.cut = PASSING_GLB_CUT + " && !triggerObjectMatchesByFilter('hltL1MuOpenL1Filtered0').empty() && " + TRACK_CUTS;

##     ___        _               _     _   _ _     _            
##    / _ \ _   _| |_ _ __  _   _| |_  | | | (_)___| |_ ___  ___ 
##   | | | | | | | __| '_ \| | | | __| | |_| | / __| __/ _ \/ __|
##   | |_| | |_| | |_| |_) | |_| | |_  |  _  | \__ \ || (_) \__ \
##    \___/ \__,_|\__| .__/ \__,_|\__| |_| |_|_|___/\__\___/|___/
##                   |_|                                         
##   
process.TFileService = cms.Service("TFileService", fileName = cms.string("tnpJPsi_Data.root"))


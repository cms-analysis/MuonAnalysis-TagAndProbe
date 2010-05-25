### Here we load both the skim stuff and the reco stuff

import FWCore.ParameterSet.Config as cms


process = cms.Process("TagProbe")

process.load('FWCore.MessageService.MessageLogger_cfi')
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.source = cms.Source("PoolSource", 
    fileNames = cms.untracked.vstring(
	#'/store/relval/CMSSW_3_5_4/RelValJpsiMM/GEN-SIM-RECO/START3X_V24-v1/0004/1244DD00-2D2C-DF11-8FA3-002618943977.root',
	#'/store/relval/CMSSW_3_5_4/RelValJpsiMM/GEN-SIM-RECO/START3X_V24-v1/0004/1217EA49-A52B-DF11-9081-00304867BFAE.root',
	#'/store/relval/CMSSW_3_5_4/RelValJpsiMM/GEN-SIM-RECO/START3X_V24-v1/0003/BE67D561-A02B-DF11-9FE7-001A928116B4.root',
	#'/store/relval/CMSSW_3_5_4/RelValJpsiMM/GEN-SIM-RECO/START3X_V24-v1/0003/B892E028-9F2B-DF11-96DD-003048679164.root',
	#'/store/relval/CMSSW_3_5_4/RelValJpsiMM/GEN-SIM-RECO/START3X_V24-v1/0003/5C692C2C-9E2B-DF11-B817-002618943869.root',
        #'root://pcmssd12.cern.ch//data/mc/7TeV/Summer09/JPsiMuMu_AODSIM_72B0D83F-908B-DE11-9A88-001D09646131.root',
        #'root://pcmssd12.cern.ch//data/mc/7TeV/Summer09/ppMuX_AODSIM_EAB08D8C-78A6-DE11-88FD-001AA009562B.root',
        'root://pcmssd12.cern.ch//data/gpetrucc/PYTHIA6_7TeV_JPsiWithFSR_MC_3XY_V16_1E31_HLT_12.root'
    ),
    duplicateCheckMode = cms.untracked.string("noDuplicateCheck")
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )    

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.load("Geometry.CommonDetUnit.globalTrackingGeometry_cfi")
process.load("TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorAny_cfi")
process.load("TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorAlong_cfi")
process.load("TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorOpposite_cfi")
process.GlobalTag.globaltag = cms.string('START36_V8::All')

process.load("MuonAnalysis.TagAndProbe.skim.skimJPsi_cff")
process.load("MuonAnalysis.TagAndProbe.jpsi.tp_from_skim_Tracking_cff")
process.load("MuonAnalysis.TagAndProbe.jpsi.tp_from_skim_MuonID_cff")
process.load("MuonAnalysis.TagAndProbe.jpsi.tp_from_skim_Trigger_cff")

process.tagAndProbe = cms.Path( 
    process.slimAOD *
    process.patMuonSequence *
    process.tnpCommonSequence    *
    ( process.tnpSequenceTracking  +
      process.tnpSequenceMuonID    +
      process.tnpSequenceTrigger   )
)

from MuonAnalysis.TagAndProbe.jpsi.tp_from_skim_common_cff import addDiMuonSeparationVariables
addDiMuonSeparationVariables(process, process.tnpSequenceMuonID, process.histoMuFromTk)

## Now, if this is unskimmed then we can turn on muon efficiencies
process.histoMuFromTk.makeMCUnbiasTree  = True
process.histoMuFromCal.makeMCUnbiasTree = True
process.histoTrigger.makeMCUnbiasTree   = True
process.histoTracking.makeMCUnbiasTree  = True
## Also, if this is Summer09 we have to change the trigger table
if False:
    from MuonAnalysis.TagAndProbe.skim.skimJPsi_cff import Summer09_Trigger
    Summer09_Trigger(process)

##     ___        _               _     _   _ _     _            
##    / _ \ _   _| |_ _ __  _   _| |_  | | | (_)___| |_ ___  ___ 
##   | | | | | | | __| '_ \| | | | __| | |_| | / __| __/ _ \/ __|
##   | |_| | |_| | |_| |_) | |_| | |_  |  _  | \__ \ || (_) \__ \
##    \___/ \__,_|\__| .__/ \__,_|\__| |_| |_|_|___/\__\___/|___/
##                   |_|                                         
##   
process.TFileService = cms.Service("TFileService", fileName = cms.string("tnpJPsi_fromAOD.root"))



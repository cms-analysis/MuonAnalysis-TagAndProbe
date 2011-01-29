import FWCore.ParameterSet.Config as cms

process = cms.Process("TagProbeTree")

process.load('FWCore.MessageService.MessageLogger_cfi')
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.MessageLogger.cerr.FwkReport.reportEvery = 100


process.source = cms.Source("PoolSource", 
    fileNames = cms.untracked.vstring(
#   'file:/tmp/botta/store/user/fat/MuOnia/Run2010B-Nov4ReReco_v1-Onia2MuMu-v6/52b0ce2a29246dfcd7fb39d814d6fb33/onia2MuMuPAT_9_1_BdS.root'
#    'file:///castor/cern.ch/user/h/hwoehri/cmst3/Fall10/JPsiToMuMu_2MuPEtaFilter_7TeV-pythia6-evtgen/PAT/onia2MuMuPAT_9_1_PfC.root'
    'rfio:/castor/cern.ch/user/h/hwoehri/cmst3/Fall10/JPsiToMuMu_2MuPEtaFilter_7TeV-pythia6-evtgen/PAT/onia2MuMuPAT_9_1_PfC.root' #MC
#    'rfio:/castor/cern.ch/user/h/hwoehri/cmst3/Data2010/4thNovReReco/PAT/PAT-MuOnia-Run2010B-Nov4ReReco-Onia2MuMu-v6-onia2MuMuPAT_9_1_BdS.root' #data
    ),
)


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1))
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000))

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
#process.GlobalTag.globaltag = cms.string('FT_R_38X_V14A::All') #ReReco4thNov
process.GlobalTag.globaltag = cms.string('START38_V12::All') #MC

##==== Merge CaloMuons and Tracks into the collection of reco::Muons ====> it is done in the previous step
##==== Make the PatMuonsWithTrigger =====> it is done in the previous step
##==== Make the tag muons, make the probe muons, make the tag&probe pairs  ====> it is done in the previous step


## Here we want to access to include all the variables we want to ask the probe to 
from MuonAnalysis.TagAndProbe.common_variables_cff import *
process.load("MuonAnalysis.TagAndProbe.common_modules_cff")

process.tagMuonsMCMatch = cms.EDProducer("MCTruthDeltaRMatcherNew",
    src = cms.InputTag("tagMuons"),
    matched = cms.InputTag("genMuons"),
    pdgId = cms.vint32(13),
    distMin = cms.double(0.3),
)

process.probeMuonsMCMatch = process.tagMuonsMCMatch.clone(src = "probeMuons")

process.tpTree = cms.EDAnalyzer("TagProbeFitTreeProducer",
                                # choice of tag and probe pairs, and arbitration
                                tagProbePairs = cms.InputTag("tpPairs"),
                                arbitration   = cms.string("OneProbe"),
                                # probe variables 
                                variables = cms.PSet(AllVariables,
                                                     dxyPVdzmin       = cms.InputTag("moreProbeInfo","dxyPVdzmin"),
                                                     ),
                               
                                flags = cms.PSet(
                                                 TrackQualityFlags,
                                                 MuonIDFlags,
                                                 HighPtTriggerFlags,
                                                 LowPtTriggerFlagsPhysics,
                                                 LowPtTriggerFlagsEfficiencies,
                                                 ## Acceptance definition
                                                 Acc_JPsi = cms.string("(abs(eta) <= 1.3 && pt > 3.3) || (1.3 < abs(eta) <= 2.2 && p > 2.9) || (2.2 < abs(eta) <= 2.4  && pt > 0.8)"),
                                                 ),
                                tagVariables = cms.PSet(
                                                   pt  = cms.string('pt'),
                                                   eta = cms.string('eta'),
                                                   nVertices = cms.InputTag("nverticesModule"),
                                                  ),
                                tagFlags = cms.PSet(LowPtTriggerFlagsPhysics,
                                                    LowPtTriggerFlagsEfficiencies,
                                                    ),
                                
                                
                                isMC             = cms.bool(True),
                                tagMatches       = cms.InputTag("tagMuonsMCMatch"),
                                probeMatches     = cms.InputTag("probeMuonsMCMatch"),
                                motherPdgId      = cms.vint32(443),
                                makeMCUnbiasTree       = cms.bool(True),
                                checkMotherInUnbiasEff = cms.bool(True),
                                allProbes              = cms.InputTag("probeMuons"),
                                
                                )



#from HLTrigger.HLTfilters.hltHighLevelDev_cfi import hltHighLevelDev
#process.fastFilter = hltHighLevelDev.clone(HLTPaths = ['HLT_Mu3_Track0_Jpsi'], HLTPathsPrescales = [1])
from HLTrigger.HLTfilters.hltHighLevel_cfi import hltHighLevel
process.fastFilter = hltHighLevel.clone(HLTPaths = ['HLT_Mu*_L2Mu0', 'HLT_Mu3_Track*_Jpsi*', 'HLT_Mu5_Track*_Jpsi*'],
                                        andOr = cms.bool(True),  # True: OR - False: AND
                                        throw = cms.bool(False), # Don't throw exception if trigger name is unkown
                                        )




process.tnpSimpleSequence = cms.Sequence(

    process.tagMuonsMCMatch    +
    process.probeMuonsMCMatch  +
    process.moreProbeInfo      +
    process.nverticesModule    +
    process.tagProbeSeparation +
    process.tpTree
)

process.tagAndProbe = cms.Path( 

    process.fastFilter         +
    process.tnpSimpleSequence             
  
)

process.TFileService = cms.Service("TFileService", fileName = cms.string("tnpTree_FromOniaMuMuPAT_Test.root"))















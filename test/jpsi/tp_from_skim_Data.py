import FWCore.ParameterSet.Config as cms

process = cms.Process("TagProbe")

process.load('FWCore.MessageService.MessageLogger_cfi')
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.source = cms.Source("PoolSource", 
    fileNames = cms.untracked.vstring(
        'root://pcmssd12.cern.ch//data/gpetrucc/7TeV/tnp/new.run139368to139400/skimJPsiLoose_Data_1_1_NDG.root'
    )
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(5000) )    

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Geometry.CommonDetUnit.globalTrackingGeometry_cfi")
process.load("TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorAny_cfi")
process.load("TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorAlong_cfi")
process.load("TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorOpposite_cfi")
process.GlobalTag.globaltag = cms.string('GR_R_36X_V12::All')

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

## Add Run and LS info in each tree entry
for K,V in allTPTreeProducers(process): 
    V.addRunLumiInfo = cms.bool(True)

from MuonAnalysis.TagAndProbe.jpsi.tp_from_skim_common_cff import addDiMuonSeparationVariables
addDiMuonSeparationVariables(process, process.tnpSequenceTrigger, process.histoTrigger)
addCountVariables(process, process.tnpSequenceTrigger,  process.histoTrigger)
addCountVariables(process, process.tnpSequenceMuonID,   process.histoMuFromTk)
addCountVariables(process, process.tnpSequenceMuonID,   process.histoMuFromCal)
addCountVariables(process, process.tnpSequenceTracking, process.histoTracking)
addCountVariables(process, process.tnpSequenceTracking, process.histoTrackingHp)

from MuonAnalysis.TagAndProbe.jpsi.tp_from_skim_Trigger_cff import ReMatchL1
ReMatchL1(process)

#process.load("UserCode.GPetrucc.edmLumi_cfi") ## see twiki UserCodeGPetrucc if you want this

##     ___        _               _     _   _ _     _            
##    / _ \ _   _| |_ _ __  _   _| |_  | | | (_)___| |_ ___  ___ 
##   | | | | | | | __| '_ \| | | | __| | |_| | / __| __/ _ \/ __|
##   | |_| | |_| | |_| |_) | |_| | |_  |  _  | \__ \ || (_) \__ \
##    \___/ \__,_|\__| .__/ \__,_|\__| |_| |_|_|___/\__\___/|___/
##                   |_|                                         
##   
process.TFileService = cms.Service("TFileService", fileName = cms.string("tnpJPsi_Data.root"))
## PVT (all detectors good):
##    GET LUMIS FROM https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions10/7TeV/
##    OR FROM /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification
## B-PAG (only tracker & muons):
##    GET FROM HeavyFlavorAnalysis/Onia2MuMu/certification/7TeV/Collisions10

## Now we use B-PAG
process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange()

# ---- CERTIFICATION (Signed off 04.07.2010) ------
# QFLAGS=Pix:GOOD,Strip:GOOD,Csc:GOOD,Dt:GOOD,Rpc:GOOD,Muon:GOOD
# DCS=Bpix,Fpix,Tibtid,TecM,TecP,Tob,Dtminus,Dtplus,Dt0,CscMinus,CscPlus,Rpc

## Jun14th Reprocessing
process.source.lumisToProcess += [
    '132440:157-132440:401',
    '132596:382-132596:453',
    '132597:1-132597:48',
    '132598:1-132598:188',
    '132599:1-132599:538',
    '132601:1-132601:207',
    '132601:209-132601:259',
    '132601:261-132601:1131',
    '132602:1-132602:83',
    '132605:1-132605:968',
    '132606:1-132606:37',
    '132656:1-132656:140',
    '132658:1-132658:177',
    '132659:1-132659:84',
    '132661:1-132661:130',
    '132662:1-132662:130',
    '132662:132-132662:165',
    '132716:220-132716:640',
    '132959:1-132959:417',
    '132960:1-132960:190',
    '132961:1-132961:427',
    '132965:1-132965:107',
    '132968:1-132968:173',
    '133029:101-133029:115',
    '133029:129-133029:350',
    '133031:1-133031:18',
    '133034:131-133034:325',
    '133035:1-133035:306',
    '133036:1-133036:225',
    #'133038:1-133038:70',  # Trigger
    '133046:1-133046:43',
    '133046:45-133046:323',
    '133082:1-133082:608',
    '133158:65-133158:786',
    #'133320:104-133320:196', # Trigger
    '133321:1-133321:383',
    '133446:105-133446:273',
    '133448:1-133448:516',
    '133450:1-133450:658',
    '133474:1-133474:95',
    '133474:157-133474:189',
    '133483:94-133483:591',
    '133483:652-133483:658',
    '133874:166-133874:875',
    '133875:1-133875:49',
    '133876:1-133876:330',
    '133877:1-133877:77',
    '133877:82-133877:231',
    '133877:236-133877:1997',
    '133881:1-133881:562',
    '133885:1-133885:728',
    '133927:1-133927:57',
    '133928:1-133928:645',
    '134721:294-134721:1340',
    '134725:231-134725:704',
    '134725:1096-134725:1122',
    '134725:1142-134725:1315',
    '135059:59-135059:67',
    '135149:294-135149:3382',
    '135175:55-135175:1082',
    '135445:168-135445:1827',
    '135521:60-135521:524',
    '135523:1-135523:109',
    '135523:113-135523:225',
    '135525:1-135525:457',
    '135528:1-135528:1436',
    '135534:83-135534:90',
    '135535:75-135535:246',
    '135573:102-135573:163',
    '135575:2-135575:638',
    '135575:645-135575:1253',
    '135735:31-135735:333',
    '136033:138-136033:1117',
    '136035:1-136035:259',
    '136066:181-136066:1184',
    '136080:249-136080:262',
    #'136082:1-136082:506',  # Trigger
    '136082:1-136082:422',
    '136082:477-136082:506',
    '136087:250-136087:315',
    '136087:335-136087:354',
    '136088:1-136088:262',
    '136097:1-136097:91',
    '136098:1-136098:25',
    '136100:1-136100:1187',
    '136119:1-136119:36',
    '136294:1-136294:60',
    '136297:1-136297:102',
    '137027:98-137027:200',
    '137028:1-137028:338',
    '137028:341-137028:484',
]

## Prompt Reco
process.source.lumisToProcess += [
    #'138562:1-138562:20', # Trigger
    #'138563:1-138563:12', # Trigger
    '138564:1-138564:19',
    #'138565:1-138565:22', # Trigger
    #'138570:1-138570:33', # Trigger
    '138571:1-138571:23',
    '138572:1-138572:238',
    '138737:1-138737:88',
    '138738:1-138738:94',
    '138739:1-138739:186',
    '138742:1-138742:54',
    '138744:1-138744:39',
    '138745:1-138745:17',
    '138746:1-138746:145',
    '138747:1-138747:131',
    '138749:1-138749:244',
    '138750:1-138750:726',
    '138751:1-138751:147',
    '138919:62-138919:168',
    '138920:1-138920:74',
    '138921:1-138921:203',
    '138923:1-138923:17',
    '138924:1-138924:86',
    '138937:1-138937:47',
    '138939:1-138939:36',
    '139020:227-139020:617',
    '139096:193-139096:280',
    '139098:1-139098:201',
    '139100:1-139100:315',
    '139103:7-139103:449',
    '139195:9-139195:92',
    '139239:164-139239:270',
    '139347:214-139347:526',
    '139356:175-139356:213',
    '139364:1-139364:89',
    '139365:1-139365:259',
    '139368:1-139368:174',
    '139370:1-139370:645',
    '139372:20-139372:202',
    '139375:1-139375:47',
    '139399:75-139399:137',
    '139400:1-139400:144',
    '139407:1-139407:1329',
    '139411:1-139411:135',
    '139457:18-139457:88',
    '139458:1-139458:410',
    '139459:1-139459:395',
]




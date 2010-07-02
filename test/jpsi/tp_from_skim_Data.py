import FWCore.ParameterSet.Config as cms

process = cms.Process("TagProbe")

process.load('FWCore.MessageService.MessageLogger_cfi')
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.source = cms.Source("PoolSource", 
    fileNames = cms.untracked.vstring(
        'root://pcmssd12.cern.ch//data/gpetrucc/7TeV/tnp/DATA_skimJPsiLoose_fromMay6th_v1.root',
        #'rfio:/castor/cern.ch/user/g/gpetrucc/7TeV/DATA/DATA_skimJPsiLoose_fromApr20MuonSkim-v2.root',
        #'rfio:/castor/cern.ch/user/g/gpetrucc/7TeV/DATA/DATA_skimJPsiLoose_fromMuonSkimV9_upToApr28-v2.root',
        #'root://pcmssd12.cern.ch//data/gpetrucc/7TeV/tnp/DATA_skimJPsiLoose_fromApr20MuonSkim-v2.root',
        #'root://pcmssd12.cern.ch//data/gpetrucc/7TeV/tnp/DATA_skimJPsiLoose_fromMuonSkimV9_upToApr28-v2.root',
        #'root://pcmssd12.cern.ch//data/gpetrucc/7TeV/tnp/DATA_skimJPsiLoose_fromPromptReco_run13509to135175.part1.root',
        #'root://pcmssd12.cern.ch//data/gpetrucc/7TeV/tnp/DATA_skimJPsiLoose_fromPromptReco_run13509to135175.part2.root',
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
process.GlobalTag.globaltag = cms.string('GR_R_36X_V11::All')

process.load("MuonAnalysis.TagAndProbe.jpsi.tp_from_skim_Tracking_cff")
process.load("MuonAnalysis.TagAndProbe.jpsi.tp_from_skim_MuonID_cff")
process.load("MuonAnalysis.TagAndProbe.jpsi.tp_from_skim_Trigger_cff")

process.tagAndProbe = cms.Path( 
    process.tnpCommonSequence    *
    ( process.tnpSequenceTracking  +
      process.tnpSequenceMuonID    +
      #process.tnpSequenceMuonIDVtx +
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
addDiMuonSeparationVariables(process, process.tnpSequenceMuonID,  process.histoMuFromTk)
addDiMuonSeparationVariables(process, process.tnpSequenceMuonID,  process.histoMuFromCal)

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

## May6th reprocessing
process.source.lumisToProcess += [
    '132440:157-132440:378',
    '132442:1-132442:83',
    '132442:86-132442:124',
    '132442:129-132442:133',
    '132442:136-132442:169',
    '132442:173-132442:214',
    '132442:218-132442:251',
    '132442:254-132442:254',
    '132442:256-132442:256',
    '132442:258-132442:258',
    '132596:382-132596:382',
    '132596:384-132596:447',
    '132597:1-132597:6',
    '132597:19-132597:40',
    '132598:1-132598:79',
    '132598:90-132598:176',
    '132599:1-132599:437',
    '132601:1-132601:207',
    '132601:209-132601:259',
    '132601:261-132601:1107',
    '132602:1-132602:70',
    '132605:1-132605:522',
    '132605:526-132605:814',
    '132605:816-132605:829',
    '132605:831-132605:867',
    '132605:896-132605:942',
    '132606:1-132606:26',
    '132656:1-132656:111',
    '132658:1-132658:51',
    '132658:56-132658:120',
    '132658:127-132658:148',
    '132659:1-132659:76',
    '132661:1-132661:116',
    '132662:1-132662:9',
    '132662:25-132662:74',
    '132716:220-132716:436',
    '132716:440-132716:487',
    '132716:491-132716:586',
    '132959:326-132959:334',
    '132960:1-132960:124',
    '132961:1-132961:222',
    '132961:226-132961:230',
    '132961:237-132961:381',
    '132965:1-132965:68',
    '132968:1-132968:67',
    '132968:75-132968:169',
    '133029:101-133029:115',
    '133029:129-133029:332',
    '133031:1-133031:18',
    '133034:132-133034:287',
    '133035:1-133035:63',
    '133035:67-133035:302',
    '133036:1-133036:222',
    '133046:1-133046:43',
    '133046:45-133046:210',
    '133046:213-133046:227',
    '133046:229-133046:323',
    '133082:1-133082:608',
    '133158:65-133158:786',
    '133321:1-133321:383',
    '133446:105-133446:266',
    '133448:1-133448:484',
    '133450:1-133450:658',
    '133474:1-133474:95',
    '133483:94-133483:591',
    '133874:166-133874:297',
    '133874:299-133874:721',
    '133874:724-133874:864',
    '133875:1-133875:37',
    '133876:1-133876:315',
    '133877:1-133877:77',
    '133877:82-133877:104',
    '133877:113-133877:231',
    '133877:236-133877:294',
    '133877:297-133877:437',
    '133877:439-133877:622',
    '133877:624-133877:853',
    '133877:857-133877:1472',
    '133877:1474-133877:1931',
    '133881:1-133881:551',
    '133885:1-133885:728',
    '133927:1-133927:44',
    '133928:1-133928:645'
]

## Prompt reco, certified, up to run 136xxx
process.source.lumisToProcess += [
    '135059:59-135059:67',
    '135149:297-135149:937',
    '135149:942-135149:993',
    '135149:995-135149:1031',
    '135149:1033-135149:1098',
    '135149:1102-135149:2269',
    '135149:2274-135149:2524',
    '135149:2528-135149:2713',
    '135149:2715-135149:3098',
    '135149:3100-135149:3102',
    '135149:3105-135149:3179',
    '135149:3182-135149:3303',
    '135149:3305-135149:3381',
    '135175:55-135175:545',
    '135175:548-135175:1046',
    '135445:168-135445:200',
    '135445:203-135445:272',
    '135445:275-135445:458',
    '135445:474-135445:832',
    '135445:834-135445:1067',
    '135445:1069-135445:1388',
    '135445:1391-135445:1629',
    '135445:1631-135445:1827',
    '135521:60-135521:108',
    '135521:110-135521:359',
    '135521:361-135521:488',
    '135523:1-135523:64',
    '135523:66-135523:109',
    '135523:113-135523:211',
    '135525:1-135525:143',
    '135525:145-135525:381',
    '135525:384-135525:435',
    '135525:437-135525:452',
    '135528:1-135528:91',
    '135528:94-135528:142',
    '135528:145-135528:308',
    '135528:310-135528:454',
    '135528:456-135528:606',
    '135528:608-135528:609',
    '135528:611-135528:770',
    '135528:773-135528:776',
    '135528:779-135528:912',
    '135528:915-135528:1082',
    '135528:1084-135528:1213',
    '135528:1215-135528:1436',
    '135535:75-135535:232',
    '135537:39-135537:69',
    '135573:102-135573:110',
    '135573:113-135573:118',
    '135573:120-135573:155',
    '135575:2-135575:241',
    '135575:243-135575:264',
    '135575:266-135575:638',
    '135575:645-135575:1161',
    '135575:1163-135575:1253',
    '135735:31-135735:42',
    '135735:44-135735:149',
    '135735:151-135735:234',
    '135735:236-135735:320',
    '136033:138-136033:1117',
    '136035:1-136035:259',
    '136066:181-136066:1184',
    '136080:249-136080:262',
    '136082:1-136082:506',
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
    '138562:1-138562:20',
    '138563:1-138563:12',
    '138564:1-138564:19',
    '138565:1-138565:22',
    '138570:1-138570:33',
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
]





import FWCore.ParameterSet.Config as cms

### USAGE:
###    cmsRun fitMuonID.py <scenario>
### scenarios:
###   - data_all (default)  
###   - signal_mc

import sys
args = sys.argv[1:]
if (sys.argv[0] == "cmsRun"): args =sys.argv[2:]
scenario = "data_all"
if len(args) > 0: scenario = args[0]
print "Will run scenario ", scenario 

process = cms.Process("TagProbe")

process.load('FWCore.MessageService.MessageLogger_cfi')

process.source = cms.Source("EmptySource")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

Template = cms.EDAnalyzer("TagProbeFitTreeAnalyzer",
   NumCPU = cms.uint32(1),
   SaveWorkspace = cms.bool(False),
   Variables = cms.PSet(
     mass = cms.vstring("Tag-muon Mass", "2.8", "3.35", "GeV/c^{2}"),
     p  = cms.vstring("muon p", "0", "1000", "GeV/c"),
     pt = cms.vstring("muon p_{T}", "0", "1000", "GeV/c"),
     eta    = cms.vstring("muon #eta", "-2.5", "2.5", ""),
     abseta = cms.vstring("muon |#eta|", "0", "2.5", ""),
     tkTrackerLay = cms.vstring("number of valid tracker layers", "-0.5", "30.5", ""),
     tkPixelLay   = cms.vstring("number of valid pixel layers", "-0.5", "12.5", ""),
     numberOfMatchedStations = cms.vstring("number of matched stations", "-0.5", "7.5", ""),
     dB = cms.vstring("dB", "-1000", "1000", ""),
     dzPV = cms.vstring("dzPV", "-1000", "1000", ""),
     tag_nVertices = cms.vstring("Number of vertices", "0", "999", ""),
     #pair_dphiVtxTimesQ = cms.vstring("q1 * (#phi1-#phi2)", "-6", "6", ""),
     pair_drM1   = cms.vstring("#Delta R(Station 1)", "-99999", "999999", "rad"),
     #pair_distM1 = cms.vstring("q", "-99999", "999999", ""),
     run = cms.vstring("Run number", "146000", "999999", ""),
     ),
                          
   Categories = cms.PSet(
     Glb   = cms.vstring("Global", "dummy[pass=1,fail=0]"),
     GlbPT = cms.vstring("Global", "dummy[pass=1,fail=0]"),
     TMA   = cms.vstring("Global", "dummy[pass=1,fail=0]"),
     Tight2012 = cms.vstring("Tight",   "dummy[pass=1,fail=0]"),
     PF        = cms.vstring("PF Muon", "dummy[pass=1,fail=0]"),
     Mu5_Track2_Jpsi_TK = cms.vstring("ProbeTrigger_Track0", "dummy[pass=1,fail=0]"),
     Mu5_Track3p5_Jpsi_TK = cms.vstring("ProbeTrigger_Track0", "dummy[pass=1,fail=0]"),
     Mu7_Track7_Jpsi_TK = cms.vstring("ProbeTrigger_Track0", "dummy[pass=1,fail=0]"),
     tag_Mu5_Track2_Jpsi_MU = cms.vstring("ProbeTrigger_Track0", "dummy[pass=1,fail=0]"),
     tag_Mu5_Track3p5_Jpsi_MU = cms.vstring("ProbeTrigger_Track0", "dummy[pass=1,fail=0]"),
     tag_Mu7_Track7_Jpsi_MU = cms.vstring("ProbeTrigger_Track0", "dummy[pass=1,fail=0]"),
     #mcTrue = cms.vstring("MC true", "dummy[true=1,false=0]"),
     #Mu5_Track0_Jpsi_TK = cms.vstring("ProbeTrigger_Track0", "dummy[pass=1,fail=0]"),
     #Mu3_Track3_Jpsi_TK = cms.vstring("ProbeTrigger_Track0", "dummy[pass=1,fail=0]"),
     #Mu5_Track5_Jpsi_TK = cms.vstring("ProbeTrigger_Track0", "dummy[pass=1,fail=0]"),
     #tag_Mu5_Track0_Jpsi_MU = cms.vstring("ProbeTrigger_Track0", "dummy[pass=1,fail=0]"),
     #tag_Mu3_Track3_Jpsi_MU = cms.vstring("ProbeTrigger_Track0", "dummy[pass=1,fail=0]"),
     #tag_Mu5_Track5_Jpsi_MU = cms.vstring("ProbeTrigger_Track0", "dummy[pass=1,fail=0]"),
     ),
                          

   Expressions = cms.PSet(
     ZZLooseVar = cms.vstring("", "Glb==1 || numberOfMatchedStations > 0", "Glb", "numberOfMatchedStations"),
     Tight2012noIPVar = cms.vstring("", "GlbPT==1 && numberOfMatchedStations > 1 && tkTrackerLay > 5 && tkPixelLay > 0", 
                    "GlbPT", "numberOfMatchedStations", "tkTrackerLay", "tkPixelLay"),
     Tight2012withIPVar = cms.vstring("", "Tight2012==1 && abs(dzPV) < 0.5", 
                    "Tight2012", "dzPV"),
   ),
                          
   Cuts = cms.PSet(
     ZZLoose         = cms.vstring("ZZLoose",  "ZZLooseVar",   "0.5"),
     Tight2012noIP   = cms.vstring("Tight2012noIP",  "Tight2012noIPVar",   "0.5"),
     Tight2012withIP = cms.vstring("Tight2012withIP","Tight2012withIPVar", "0.5"),
   ),
                          
                          
   PDFs = cms.PSet(
     gaussPlusExpo = cms.vstring(
     #"CBShape::signal(mass, mean[3.1,3.0,3.2], sigma[0.05,0.02,0.06], alpha[3., 0.5, 5.], n[1, 0.1, 100.])",
     #"Chebychev::backgroundPass(mass, {cPass[0,-0.5,0.5], cPass2[0,-0.5,0.5]})",
     #"Chebychev::backgroundFail(mass, {cFail[0,-0.5,0.5], cFail2[0,-0.5,0.5]})",
     "Gaussian::signal(mass, mean[3.1,3.0,3.2], sigma[0.05,0.02,0.1])",
     "Exponential::backgroundPass(mass, lp[0,-5,5])",
     "Exponential::backgroundFail(mass, lf[0,-5,5])",
     "efficiency[0.9,0,1]",
     "signalFractionInPassing[0.9]"
     )
   ),
                          
   binnedFit = cms.bool(True),
   binsForFit = cms.uint32(40),
   saveDistributionsPlot = cms.bool(False),
                          
   Efficiencies = cms.PSet(), # will be filled later
)

# pick muons that bend apart from each other
SEPARATED = cms.PSet(pair_drM1 = cms.vdouble(0.5,10))

PT_ETA_BINS = cms.PSet(
   SEPARATED,
   pt     = cms.vdouble(  2, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.5, 5.0, 5.5, 6.0, 7.5, 9.0, 11.0, 15.0, 20.0), # --> Gio's thesis bins
   abseta = cms.vdouble(  0.0, 1.2, 2.4),
   p      = cms.vdouble(  3.5, 999),
)

PT_ETA_BINS2 = cms.PSet(SEPARATED,
   pt     = cms.vdouble(  3.0, 5.0, 7.5, 10.0, 15.0, 20.0), # --> SF bins
   abseta = cms.vdouble(  0.0, 1.2, 2.4),
   p      = cms.vdouble(  3.5, 999),
)

VTX_BINS = cms.PSet(SEPARATED,
                    abseta = cms.vdouble(0.0, 2.4),
                    pt     = cms.vdouble(7.0, 20.0),
                    tag_nVertices = cms.vdouble(0.5,2.5,4.5,6.5,8.5,10.5,12.5,16.5,20.5),
                    )


process.TnP_MuonID = Template.clone(
    InputFileNames = cms.vstring(),
    InputTreeName = cms.string("fitter_tree"),
    InputDirectoryName = cms.string("tpTree"),
    OutputFileName = cms.string("TnPJPsi_MuonID_%s.root" % scenario),
    Efficiencies = cms.PSet(),
)
if "_weight" in scenario:
    process.TnP_MuonID.WeightVariable = cms.string("weight")
    process.TnP_MuonID.Variables.weight = cms.vstring("weight","0","10","")


if "data2011" in scenario:
    process.TnP_MuonID.InputFileNames = [ 
        "root://eoscms//eos/cms/store/cmst3/user/botta/TnPtrees/tnpJPsi_Run2011Amay10.root",
        "root://eoscms//eos/cms/store/cmst3/user/botta/TnPtrees/tnpJPsi_Run2011Aaug05.root",
        "root://eoscms//eos/cms/store/cmst3/user/botta/TnPtrees/tnpJPsi_Run2011Av6.root",
        "root://eoscms//eos/cms/store/cmst3/user/botta/TnPtrees/tnpJPsi_Run2011Bv1.root",
    ]
if "mc2011" in scenario:
    process.TnP_MuonID.InputFileNames = [ 
        "root://eoscms//eos/cms/store/cmst3/user/botta/TnPtrees/tnpJPsi_JPsiToMuMu_Summer11.root",
    ]
if "data2012" in scenario:
    process.TnP_MuonID.InputFileNames = [ 
        "root://eoscms//eos/cms/store/cmst3/user/botta/TnPtrees/tnpJPsi_Data.190456-193557.root",
        "root://eoscms//eos/cms/store/cmst3/user/botta/TnPtrees/tnpJPsi_Data_193557-193752.root",
        "root://eoscms//eos/cms/store/cmst3/user/botta/TnPtrees/tnpJPsi_Data_194108-194479.root",
        "root://pcmssd12//data/gpetrucc/8TeV/tnp/All/tnpJPsi_Data_194480-195016.root",
        "root://pcmssd12//data/gpetrucc/8TeV/tnp/All/tnpJPsi_Data_195017-195775.root",
        "root://pcmssd12//data/gpetrucc/8TeV/tnp/All/tnpJPsi_Data_195776-196531.root",
    ]
if "mc2012" in scenario:
    process.TnP_MuonID.InputFileNames = [ 
        #"root://pcmssd12//data/gpetrucc/8TeV/tnp/All/tnpJPsi_BsToJPsiPhi_Summer12_PU_S7_START52_V5_GoodEvents.root",
        #"root://pcmssd12//data/gpetrucc/8TeV/tnp/All/tnpJPsi_BuToJPsiK_Summer12_PU_S7_START52_V9_GoodEvents.root",
        "root://pcmssd12//data/gpetrucc/8TeV/tnp/All/tnpJPsi_MC_BsToJPsiPhi_2K2MuPtEtaFilter_Summer12_PU_S7_START52_V9_v1.root",
        "root://pcmssd12//data/gpetrucc/8TeV/tnp/All/tnpJPsi_MC_BsToJPsiPhi_2K2MuFilter_Summer12_PU_S7_START52_V9_v2.root",
    ]



IDS = ["PF", "Tight2012noIP", "Tight2012withIP", "ZZLoose", "Glb", "TMA" ]
IDSABOVE = [ "Tight2012noIP", "Tight2012withIP", "ZZLoose" ]
TRIGS = [ (2,'Mu5_Track2'), (3.75,'Mu5_Track3p5'), (7.5,'Mu7_Track7') ]
if "2011" in scenario:
     TRIGS = [ (2,'Mu5_Track2'), (7.5,'Mu7_Track7') ]

ALLBINS =  [("pt_abseta",PT_ETA_BINS)]
ALLBINS +=  [("pt_abseta2",PT_ETA_BINS2)]

for ID in IDS:
     if len(args) > 1 and args[1] in IDS and ID != args[1]: continue
     for X,B in ALLBINS:
          if len(args) > 2 and X not in args[2:]: continue
          module = process.TnP_MuonID.clone(OutputFileName = cms.string("TnP_JPsi_MuonID_%s_%s_%s.root" % (scenario, ID, X)))
          for PTMIN, TRIG in TRIGS: 
               TRIGLABEL=""
               if "pt_" in X:
                    TRIGLABEL="_"+TRIG
               else:
                    if TRIG != "Mu7_Track7": continue # use only one trigger except for turn-on
               DEN = B.clone()
               if   "2011A" in scenario: DEN.run = cms.vdouble(160404,173692)
               elif "2011B" in scenario: DEN.run = cms.vdouble(175832,177791)
               if hasattr(DEN, "pt"):
                    DEN.pt = cms.vdouble(*[i for i in B.pt if i >= PTMIN])
                    if len(DEN.pt) == 1: DEN.pt = cms.vdouble(PTMIN, DEN.pt[0])
               setattr(DEN, "tag_%s_Jpsi_MU" % TRIG, cms.vstring("pass"))
               setattr(DEN,     "%s_Jpsi_TK" % TRIG, cms.vstring("pass"))
               SUCC = "pass"
               if ID in IDSABOVE: SUCC = "above"
               setattr(module.Efficiencies, ID+"_"+X+TRIGLABEL, cms.PSet(
                         EfficiencyCategoryAndState = cms.vstring(ID, SUCC),
                         UnbinnedVariables = cms.vstring("mass"),
                         BinnedVariables = DEN,
                         BinToPDFmap = cms.vstring("gaussPlusExpo")
                         ))
               #if "plateau" in X: module.SaveWorkspace = True
               ## mc efficiency, if scenario is mc 
               if "mc" in scenario:
                    setattr(module.Efficiencies, ID+"_"+X+TRIGLABEL+"_mcTrue", cms.PSet(
                              EfficiencyCategoryAndState = cms.vstring(ID, SUCC),
                              UnbinnedVariables = cms.vstring("mass"),
                              BinnedVariables = DEN.clone(mcTrue = cms.vstring("true"))
                              ))
          setattr(process, "TnP_MuonID_"+ID+"_"+X, module)        
          setattr(process, "run_"+ID+"_"+X, cms.Path(module))


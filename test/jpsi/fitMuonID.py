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
        eta = cms.vstring("muon #eta", "-2.5", "2.5", ""),
        abseta = cms.vstring("muon |#eta|", "0", "2.5", ""),
        tag_pt = cms.vstring("Tag p_{T}",    "0", "1000", "GeV/c"),
        tag_Nvertices = cms.vstring("Number of vertices", "0", "999", ""),
        tag_nVerticesDA = cms.vstring("Number of vertices", "0", "999", ""),
        pair_dphiVtxTimesQ = cms.vstring("q1 * (#phi1-#phi2)", "-6", "6", ""),
        pair_drM1   = cms.vstring("#Delta R(Station 1)", "-99999", "999999", "rad"),
        pair_distM1 = cms.vstring("q", "-99999", "999999", ""),
    ),

    Categories = cms.PSet(
        Glb   = cms.vstring("Global", "dummy[pass=1,fail=0]"),
        VBTF  = cms.vstring("VBTFLike", "dummy[pass=1,fail=0]"),
        TMOST = cms.vstring("TMOneStationTight", "dummy[pass=1,fail=0]"),
        PF    = cms.vstring("PF Muon", "dummy[pass=1,fail=0]"),
        Mu5_Track2_Jpsi_TK = cms.vstring("ProbeTrigger_Track0", "dummy[pass=1,fail=0]"),
        Mu7_Track7_Jpsi_TK = cms.vstring("ProbeTrigger_Track0", "dummy[pass=1,fail=0]"),
        tag_Mu5_Track2_Jpsi_MU = cms.vstring("ProbeTrigger_Track0", "dummy[pass=1,fail=0]"),
        tag_Mu7_Track7_Jpsi_MU = cms.vstring("ProbeTrigger_Track0", "dummy[pass=1,fail=0]"),
        mcTrue = cms.vstring("MC true", "dummy[true=1,false=0]"),
        Mu5_Track0_Jpsi_TK = cms.vstring("ProbeTrigger_Track0", "dummy[pass=1,fail=0]"),
        Mu3_Track3_Jpsi_TK = cms.vstring("ProbeTrigger_Track0", "dummy[pass=1,fail=0]"),
        Mu5_Track5_Jpsi_TK = cms.vstring("ProbeTrigger_Track0", "dummy[pass=1,fail=0]"),
        tag_Mu5_Track0_Jpsi_MU = cms.vstring("ProbeTrigger_Track0", "dummy[pass=1,fail=0]"),
        tag_Mu3_Track3_Jpsi_MU = cms.vstring("ProbeTrigger_Track0", "dummy[pass=1,fail=0]"),
        tag_Mu5_Track5_Jpsi_MU = cms.vstring("ProbeTrigger_Track0", "dummy[pass=1,fail=0]"),
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

    Efficiencies = cms.PSet(), # will be filled later
)

# pick muons that bend apart from each other
SEPARATED = cms.PSet(pair_drM1 = cms.vdouble(0.5,10))

PT_ETA_BINS = cms.PSet(SEPARATED,
    pt     = cms.vdouble(  2, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.5, 5.0, 5.5, 6.0, 7.0, 9.0, 11.0, 14.0, 17.0, 20.0),
    abseta = cms.vdouble(  0.0, 1.2, 2.4)
)

VTX_BINS = cms.PSet(SEPARATED,
    abseta = cms.vdouble(0.0, 1.2, 2.4),
    pt     = cms.vdouble(7.0, 20.0),
    tag_nVerticesDA = cms.vdouble(0.5,2.5,4.5,6.5,8.5,10.5,12.5,14.5,16.5,18.5,20.5)
)
PLATEAU_ABSETA = cms.PSet(SEPARATED,
    abseta = cms.vdouble(0.0,  1.2, 2.4),
    pt     = cms.vdouble(7.0, 30.0),
)

PT_ABSETA_WIDE = cms.PSet(SEPARATED,
    abseta = cms.vdouble(0.0, 1.2, 2.4),
    pt     = cms.vdouble(5.0, 7.0, 20.0),
)

OLDPREFIX="/data/gpetrucc/7TeV/tnp/2011.02.17/"
PREFIX="/data/gpetrucc/7TeV/tnp/2011-04/"
process.TnP_MuonID = Template.clone(
    InputFileNames = cms.vstring(
        PREFIX+'tnpJPsi_2011A_v1_MuonPhys.root',
        PREFIX+'tnpJPsi_2011A_v2_MuonPhys.root',
    ),
    InputTreeName = cms.string("fitter_tree"),
    InputDirectoryName = cms.string("tpTree"),
    OutputFileName = cms.string("TnP_MuonID_%s.root" % scenario),
    Efficiencies = cms.PSet(),
)

IDS = [ "Glb", "TMOST", "VBTF", "PF" ]
TRIGS = [ (2,'Mu5_Track2'), (7,'Mu7_Track7') ]

if "mc_39X" in scenario:
     TRIGS = [ (0,'Mu5_Track0'), (3,'Mu3_Track3'), (5,'Mu3_Track5') ]
     process.TnP_MuonID.InputFileNames = [ OLDPREFIX+"tnpJPsi_MC_Prompt_39X.crop.root" ]

ALLBINS =  [("pt_abseta",PT_ETA_BINS), ("vtx",VTX_BINS), ("plateau",PLATEAU_ABSETA)]
ALLBINS += [("pt_abseta_wide",PT_ABSETA_WIDE)]

for ID in IDS:
    if len(args) > 1 and args[1] in IDS and ID != args[1]: continue
    for X,B in ALLBINS:
        if len(args) > 2 and X not in args[2:]: continue
        module = process.TnP_MuonID.clone(OutputFileName = cms.string("TnP_MuonID_%s_%s_%s.root" % (scenario, ID, X)))
        for PTMIN, TRIG in TRIGS: 
            TRIGLABEL=""
            if "pt_" in X:
                TRIGLABEL="_"+TRIG
            else:
                if TRIG != "Mu7_Track7": continue # use only one trigger except for turn-on
            DEN = B.clone()
            if hasattr(DEN, "pt"):
                DEN.pt = cms.vdouble(*[i for i in B.pt if i >= PTMIN])
                if len(DEN.pt) == 1: DEN.pt = cms.vdouble(PTMIN, DEN.pt[0])
            setattr(DEN, "tag_%s_Jpsi_MU" % TRIG, cms.vstring("pass"))
            setattr(DEN,     "%s_Jpsi_TK" % TRIG, cms.vstring("pass"))
            #if "calomu" in scenario: DEN.Calo = cms.vstring("pass")
            setattr(module.Efficiencies, ID+"_"+X+TRIGLABEL, cms.PSet(
                EfficiencyCategoryAndState = cms.vstring(ID,"pass"),
                UnbinnedVariables = cms.vstring("mass"),
                BinnedVariables = DEN,
                BinToPDFmap = cms.vstring("gaussPlusExpo")
            ))
            i#if "plateau" in X: module.SaveWorkspace = True
            ## mc efficiency, if scenario is mc 
            if "mc" in scenario:
                setattr(module.Efficiencies, ID+"_"+X+TRIGLABEL+"_mcTrue", cms.PSet(
                    EfficiencyCategoryAndState = cms.vstring(ID,"pass"),
                    UnbinnedVariables = cms.vstring("mass"),
                    BinnedVariables = DEN.clone(mcTrue = cms.vstring("true"))
                ))
        setattr(process, "TnP_MuonID_"+ID+"_"+X, module)        
        setattr(process, "run_"+ID+"_"+X, cms.Path(module))


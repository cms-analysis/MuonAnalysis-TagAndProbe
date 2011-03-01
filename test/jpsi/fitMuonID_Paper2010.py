import FWCore.ParameterSet.Config as cms

### USAGE:
###    cmsRun fitMuonID_Paper2010.py <scenario>
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
        pair_Nvertices = cms.vstring("Number of vertices", "0", "999", ""),
        pair_dphiVtxTimesQ = cms.vstring("q1 * (#phi1-#phi2)", "-6", "6", ""),
        pair_distM1  = cms.vstring("q1 * (#phi1-#phi2)", "-99999", "999999", ""),
        tag_nVertices = cms.vstring("Number of vertices", "0", "999", ""),
    ),

    Categories = cms.PSet(
        Calo  = cms.vstring("Global", "dummy[pass=1,fail=0]"),
        Glb   = cms.vstring("Global", "dummy[pass=1,fail=0]"),
        VBTF  = cms.vstring("VBTFLike", "dummy[pass=1,fail=0]"),
        TMOST = cms.vstring("TMOneStationTight", "dummy[pass=1,fail=0]"),
        PF    = cms.vstring("PF Muon", "dummy[pass=1,fail=0]"),
        Mu3_Track0_Jpsi_TK = cms.vstring("ProbeTrigger_Track0", "dummy[pass=1,fail=0]"),
        Mu5_Track0_Jpsi_TK = cms.vstring("ProbeTrigger_Track0", "dummy[pass=1,fail=0]"),
        Mu3_Track3_Jpsi_TK = cms.vstring("ProbeTrigger_Track0", "dummy[pass=1,fail=0]"),
        Mu3_Track5_Jpsi_TK = cms.vstring("ProbeTrigger_Track0", "dummy[pass=1,fail=0]"),
        tag_Mu3_Track0_Jpsi_MU = cms.vstring("ProbeTrigger_Track0", "dummy[pass=1,fail=0]"),
        tag_Mu5_Track0_Jpsi_MU = cms.vstring("ProbeTrigger_Track0", "dummy[pass=1,fail=0]"),
        tag_Mu3_Track3_Jpsi_MU = cms.vstring("ProbeTrigger_Track0", "dummy[pass=1,fail=0]"),
        tag_Mu3_Track5_Jpsi_MU = cms.vstring("ProbeTrigger_Track0", "dummy[pass=1,fail=0]"),
        mcTrue = cms.vstring("MC true", "dummy[true=1,false=0]"),
    ),

    PDFs = cms.PSet(
        cbPlusPoly = cms.vstring(
            "CBShape::signal(mass, mean[3.1,3.0,3.2], sigma[0.05,0.02,0.06], alpha[3., 0.5, 5.], n[1, 0.1, 100.])",
            "Chebychev::backgroundPass(mass, {cPass[0,-0.5,0.5], cPass2[0,-0.5,0.5]})",
            "Chebychev::backgroundFail(mass, {cFail[0,-0.5,0.5], cFail2[0,-0.5,0.5]})",
            "efficiency[0.9,0,1]",
            "signalFractionInPassing[0.9]"
        )
    ),

    binnedFit = cms.bool(True),
    binsForFit = cms.uint32(40),

    Efficiencies = cms.PSet(), # will be filled later
)

# pick muons that bend apart from each other
SEAGULL   = cms.PSet(pair_dphiVtxTimesQ = cms.vdouble(-2,0))
SEPARATED = cms.PSet(pair_distM1 = cms.vdouble(200,1000))

PT_ETA_BINS = cms.PSet(SEAGULL, SEPARATED,
    pt     = cms.vdouble(  0.5, 1.0, 1.5, 2, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.5, 5.0, 5.5, 6.0, 7.0, 9.0, 11.0, 14.0, 17.0, 20.0),
    abseta = cms.vdouble(  0.0, 1.2, 2.4)
)
P_ETA_BINS = cms.PSet(SEAGULL, SEPARATED,
    p      = cms.vdouble(  2, 2.5, 3, 3.5, 4, 4.5, 5, 5, 6, 8, 10, 14, 20 ),
    abseta = cms.vdouble(  1.2, 2.4)
)
ETA_BINS = cms.PSet(SEAGULL, SEPARATED,
    pt  = cms.vdouble(5,100),
    eta = cms.vdouble(-2.4, -2.1, -1.6, -1.1, -0.6, 0, 0.6, 1.1, 1.6, 2.1, 2.4),
)

VTX_BINS_BARREL = cms.PSet(
    SEAGULL, SEPARATED,
    abseta = cms.vdouble(0.0, 0.9),
    pt     = cms.vdouble(6.0, 20.0),
    tag_nVertices = cms.vdouble(0.5,1.5,2.5,3.5,4.5,5.5,6.5)
)
VTX_BINS_ENDCAPS = cms.PSet(
    SEAGULL, SEPARATED,
    abseta = cms.vdouble(1.2, 2.4),
    pt     = cms.vdouble(4.0, 20.0),
    tag_nVertices = cms.vdouble(0.5,1.5,2.5,3.5,4.5,5.5,6.5)
)


PREFIX="root://lxcms06//data2/b/botta/TnP_Paper2010_Tree/"
PREFIX="/data/gpetrucc/7TeV/tnp/2011.02.17/"
process.TnP_MuonID = Template.clone(
    InputFileNames = cms.vstring(
        PREFIX+'tnpJPsi_Data_Nov4A.root',
        PREFIX+'tnpJPsi_Data_Nov4B.root',
    ),
    InputTreeName = cms.string("fitter_tree"),
    InputDirectoryName = cms.string("tpTree"),
    OutputFileName = cms.string("TnP_Paper2010_MuonID_%s.root" % scenario),
    Efficiencies = cms.PSet(),
)

if scenario.find("signal_mc") != -1:
    process.TnP_MuonID.InputFileNames = [ PREFIX+"tnpJPsi_MC_Prompt.root" ]
if scenario.find("some_mc") != -1:
    process.TnP_MuonID.InputFileNames = [ PREFIX+"tnpJPsi_MC_Prompt.crop.root" ]
if scenario.find("beauty_mc") != -1:
    process.TnP_MuonID.InputFileNames = [ PREFIX+"tnpJPsi_MC_Bp.root" ]

#IDS = [ "Glb", "TMOST", "VBTF", "PF" ]
IDS = [ "TMOST", "VBTF", "PF" ]
TRIGS = [ (0,'Mu3_Track0'), (0,'Mu5_Track0'), (3,'Mu3_Track3'), (5,'Mu3_Track5') ]
if scenario.find("calomu") != -1:
    TRIGS = [ (0,'Mu3_Track0'), (0,'Mu5_Track0') ]

ALLBINS=[("pt_abseta",PT_ETA_BINS), ("vtx_barrel",VTX_BINS_BARREL), ("vtx_endcaps",VTX_BINS_ENDCAPS)]
ALLBINS += [ ("p_abseta",P_ETA_BINS), ("pt5_eta",ETA_BINS) ]
    
for ID in IDS:
    if len(args) > 1 and ID != args[1]: continue
    module = process.TnP_MuonID.clone(OutputFileName = cms.string("TnP_Paper2010_MuonID_%s_%s.root" % (scenario, ID)))
    for X,B in ALLBINS:
        for PTMIN, TRIG in TRIGS: 
            if X.find("vtx") != -1 and TRIG != "Mu3_Track3": continue # use only one trigger for pileup studies
            if (TRIG == "Mu3_Track5" or TRIG == "Mu3_Track0") and scenario.find("mc") != -1: continue # skip trigger not in MC
            if X.find("p_abseta") != -1 or X.find("p_2d") != -1:
                DEN=B.clone()
            else:
                DEN=B.clone(pt = cms.vdouble(*[i for i in B.pt if i >= PTMIN]))
                if len(DEN.pt) == 1: DEN.pt = cms.vdouble(PTMIN, DEN.pt[0])
            setattr(DEN, "tag_%s_Jpsi_MU" % TRIG, cms.vstring("pass"))
            setattr(DEN,     "%s_Jpsi_TK" % TRIG, cms.vstring("pass"))
            if scenario.find("calomu") != -1: DEN.Calo = cms.vstring("pass")
            ## fit efficiency, unless scenario vetoes it
            if scenario.find("nofit") == -1:
                setattr(module.Efficiencies, ID+"_"+X+"_"+TRIG, cms.PSet(
                    EfficiencyCategoryAndState = cms.vstring(ID,"pass"),
                    UnbinnedVariables = cms.vstring("mass"),
                    BinnedVariables = DEN,
                    BinToPDFmap = cms.vstring("cbPlusPoly")
                ))
            ## mc efficiency, if scenario is mc 
            if scenario.find("mc") != -1:
                setattr(module.Efficiencies, ID+"_"+X+"_"+TRIG+"_mcTrue", cms.PSet(
                    EfficiencyCategoryAndState = cms.vstring(ID,"pass"),
                    UnbinnedVariables = cms.vstring("mass"),
                    BinnedVariables = DEN.clone(mcTrue = cms.vstring("true"))
                ))
    setattr(process, "TnP_MuonID_"+ID, module)        
    setattr(process, "run_"+ID, cms.Path(module))


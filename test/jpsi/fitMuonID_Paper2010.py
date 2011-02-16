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
        mass = cms.vstring("Tag-Probe Mass", "2.8", "3.35", "GeV/c^{2}"),
        pt = cms.vstring("Probe p_{T}", "0", "1000", "GeV/c"),
        abseta = cms.vstring("Probe |#eta|", "0", "2.5", ""),
        tag_pt = cms.vstring("Tag p_{T}",    "0", "1000", "GeV/c"),
        pair_Nvertices = cms.vstring("Number of vertices", "0", "999", ""),
        pair_dphiVtxTimesQ = cms.vstring("q1 * (#phi1-#phi2)", "-6", "6", ""),
        pair_distM1  = cms.vstring("q1 * (#phi1-#phi2)", "-99999", "999999", ""),
        tag_nVertices = cms.vstring("Number of vertices", "0", "999", ""),
    ),

    Categories = cms.PSet(
        Glb   = cms.vstring("Global", "dummy[pass=1,fail=0]"),
        VBTF  = cms.vstring("VBTFLike", "dummy[pass=1,fail=0]"),
        TMOST = cms.vstring("TMOneStationTight", "dummy[pass=1,fail=0]"),
        PF    = cms.vstring("PF Muon", "dummy[pass=1,fail=0]"),
        Mu5_Track0_Jpsi_TK = cms.vstring("ProbeTrigger_Track0", "dummy[pass=1,fail=0]"),
        Mu3_Track3_Jpsi_TK = cms.vstring("ProbeTrigger_Track0", "dummy[pass=1,fail=0]"),
        Mu3_Track5_Jpsi_TK = cms.vstring("ProbeTrigger_Track0", "dummy[pass=1,fail=0]"),
        tag_Mu5_Track0_Jpsi_MU = cms.vstring("ProbeTrigger_Track0", "dummy[pass=1,fail=0]"),
        tag_Mu3_Track3_Jpsi_MU = cms.vstring("ProbeTrigger_Track0", "dummy[pass=1,fail=0]"),
        tag_Mu3_Track5_Jpsi_MU = cms.vstring("ProbeTrigger_Track0", "dummy[pass=1,fail=0]"),
        mcTrue = cms.vstring("MC true", "dummy[true=1,false=0]"),
    ),

    PDFs = cms.PSet(
        cbPlusPoly = cms.vstring(
            "CBShape::signal(mass, mean[3.1,3.0,3.2], sigma[0.05,0.02,0.06], alpha[3., 0.5, 5.], n[1, 0., 100.])",
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

PT_ETA_BINS = cms.PSet(
    pt     = cms.vdouble(  0.5, 1.0, 1.5, 2, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 7.0, 9.0, 11.0, 14.0, 17.0, 20.0),
    abseta = cms.vdouble(  0.0, 1.2, 2.4)
)
PT_ETA_BINS_SEP = cms.PSet(PT_ETA_BINS, SEAGULL, SEPARATED)

VTX_BINS_BARREL = cms.PSet(
    SEAGULL, SEPARATED,
    abseta = cms.vdouble(0.0, 1.2),
    pt     = cms.vdouble(6.0, 20.0),
    tag_nVertices = cms.vdouble(0.5,1.5,2.5,3.5,4.5,5.5,6.5)
)
VTX_BINS_ENDCAPS = cms.PSet(
    SEAGULL, SEPARATED,
    abseta = cms.vdouble(1.2, 2.4),
    pt     = cms.vdouble(4.0, 20.0),
    tag_nVertices = cms.vdouble(0.5,1.5,2.5,3.5,4.5,5.5,6.5)
)


PREFIX="/data/gpetrucc/7TeV/tnp/2011.02.11/"
process.TnP_MuonID = Template.clone(
    InputFileNames = cms.vstring(PREFIX+'tnpJPsi_MuOnia_Run2010B_Nov4.root'),
    InputTreeName = cms.string("fitter_tree"),
    InputDirectoryName = cms.string("tpTree"),
    OutputFileName = cms.string("TnP_Paper2010_MuonID_%s.root" % scenario),
    Efficiencies = cms.PSet(),
)

if scenario == "signal_mc":
    process.TnP_MuonID.InputFileNames = [ PREFIX+"tnpJPsi_MC_Prompt.root" ]
if scenario == "some_mc":
    process.TnP_MuonID.InputFileNames = [ PREFIX+"tnpJPsi_MC_Prompt.crop.root" ]
if scenario == "beauty_mc":
    process.TnP_MuonID.InputFileNames = [ PREFIX+"tnpJPsi_MC_Bp.root" ]

IDS = [ "Glb", "TMOST", "VBTF", "PF" ]
TRIGS = [ (0,'Mu5_Track0'), (3,'Mu3_Track3'), (5,'Mu3_Track5') ]
ALLBINS=[("pt_abseta",PT_ETA_BINS_SEP), ("vtx_barrel",VTX_BINS_BARREL), ("vtx_endcaps",VTX_BINS_ENDCAPS)]

for ID in IDS:
    module = process.TnP_MuonID.clone(OutputFileName = cms.string("TnP_Paper2010_MuonID_%s_%s.root" % (scenario, ID)))
    for PTMIN, TRIG in TRIGS: 
        for X,B in ALLBINS:
            DEN=B.clone(pt = cms.vdouble(*[i for i in B.pt if i >= PTMIN]))
            if len(DEN.pt) == 1: DEN.pt = cms.vdouble(PTMIN, DEN.pt[0])
            setattr(DEN, "tag_%s_Jpsi_MU" % TRIG, cms.vstring("pass"))
            setattr(DEN,     "%s_Jpsi_TK" % TRIG, cms.vstring("pass"))
            setattr(module.Efficiencies, ID+"_"+X+"_"+TRIG, cms.PSet(
                EfficiencyCategoryAndState = cms.vstring(ID,"pass"),
                UnbinnedVariables = cms.vstring("mass"),
                BinnedVariables = DEN,
                BinToPDFmap = cms.vstring("cbPlusPoly")
            ))
            if scenario == "datalike_mc" or scenario == "signal_mc":
                setattr(module.Efficiencies, ID+"_"+X+"_"+TRIG+"_mcTrue", cms.PSet(
                    EfficiencyCategoryAndState = cms.vstring(ID,"pass"),
                    UnbinnedVariables = cms.vstring("mass"),
                    BinnedVariables = DEN.clone(mcTrue = cms.vstring("true"))
                ))
    setattr(process, "TnP_MuonID_"+ID, module)        
    setattr(process, "run_"+ID, cms.Path(module))


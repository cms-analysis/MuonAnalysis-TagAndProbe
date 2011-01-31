import FWCore.ParameterSet.Config as cms

import sys
args = sys.argv[1:]
if (sys.argv[0] == "cmsRun"): args =sys.argv[2:]
scenario = "data_all"
if len(args) > 0: scenario = args[0]
print "Will run scenario ", scenario 
noEta = True # if scenario == "data_all" else True
ichep = False

CONSTRAINTS = cms.PSet(
    outerValidHits = cms.vstring("pass"),
    MuX_L2Mu0_L2     = cms.vstring("pass"),
    tag_Mu5_L2Mu0_Mu = cms.vstring("pass"),
)
ONE_BIN = cms.PSet(CONSTRAINTS,
    pt = cms.vdouble( 0, 20 ),
    #abseta = cms.vdouble(0, 2.4),
    eta = cms.vdouble(-2.4, 2.4),
)
ETA_BINS = ONE_BIN.clone(
   eta = cms.vdouble(-2.4, -2.1, -1.8, -1.6, -1.4, -1.2, -1.0, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.5, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.1, 2.4)
)
ICHEP_BINS = cms.PSet(CONSTRAINTS,
    pt = cms.vdouble( 0, 20 ),
    abseta = cms.vdouble(0, 1.1, 1.6, 2.1, 2.4)
)

ETA_PHI_BINS =  ONE_BIN.clone(
   eta = cms.vdouble(*[0.4*i for i in range(-6,7)]),
   phi = cms.vdouble(*[3.1416*i/3.0 for i in range(-3,3)]), 
)

VTX_BINS = ONE_BIN.clone(
    tag_nVertices = cms.vdouble(0.5,1.5,2.5,3.5,4.5,5.5,6.5)
)

process = cms.Process("TagProbe")

process.load('FWCore.MessageService.MessageLogger_cfi')

process.source = cms.Source("EmptySource")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

Template = cms.EDAnalyzer("TagProbeFitTreeAnalyzer",
    NumCPU = cms.uint32(1),
    SaveWorkspace = cms.bool(False),

    Variables = cms.PSet(
        mass   = cms.vstring("Tag-Probe Mass", "2.0", "4.3", "GeV/c^{2}"),
        pt     = cms.vstring("Probe p_{T}", "0", "1000", "GeV/c"),
        p      = cms.vstring("Probe p", "0", "1000", "GeV/c"),
        eta    = cms.vstring("Probe #eta", "-2.5", "2.5", ""),
        abseta = cms.vstring("Probe |#eta|", "0", "2.5", ""),
        phi    = cms.vstring("Probe #phi", "-3.1416", "3.1416", ""),
        tag_pt = cms.vstring("Tag p_{T}", "2.6", "1000", "GeV/c"),
        tk_deltaR   = cms.vstring("Match #Delta R",    "0", "1000", ""),
        tk_deltaEta = cms.vstring("Match #Delta #eta", "0", "1000", ""),
        tk_deltaR_NoJPsi   = cms.vstring("Unmatch #Delta R",    "0", "1000", ""),
        tk_deltaEta_NoJPsi = cms.vstring("Unmatch #Delta #eta", "0", "1000", ""),
        tk_deltaR_NoBestJPsi   = cms.vstring("Unmatch #Delta R",    "0", "1000", ""),
        tk_deltaEta_NoBestJPsi = cms.vstring("Unmatch #Delta #eta", "0", "1000", ""),
        tag_nVertices = cms.vstring("Number of vertices", "0", "999", ""),
        run    = cms.vstring("Run Number", "132440", "999999", ""),
    ),

    Categories = cms.PSet(
        outerValidHits = cms.vstring("hasValidHits",  "dummy[pass=1,fail=0]"),
        MuX_L2Mu0_L2     = cms.vstring("probeL2Mu0",  "dummy[pass=1,fail=0]"),
        tag_Mu5_L2Mu0_Mu = cms.vstring("tagMu5L2Mu0", "dummy[pass=1,fail=0]"),
        mcTrue = cms.vstring("MC true", "dummy[true=1,false=0]"),
    ),

    Expressions = cms.PSet(),
    Cuts = cms.PSet(),

    PDFs = cms.PSet(
        gaussPlusCubic = cms.vstring(
            "Gaussian::signal(mass, mean[3.1,3.0,3.2], sigma[0.15,0.05,0.25])",
            "Chebychev::backgroundPass(mass, {c1p[0,-1,1], c2p[0,-1,1], c3p[0,-1,1]})",
            "Chebychev::backgroundFail(mass, {c1f[0,-1,1], c2f[0,-1,1], c3f[0,-1,1]})",
            "efficiency[0.9,0,1]",
            "signalFractionInPassing[0.5]"
        ),
        gaussPlusFloatCubic = cms.vstring(
            "Gaussian::signal(mass, mean[3.1,3.0,3.2], sigma[0.15,0.05,0.25])",
            "Chebychev::backgroundPass(mass, {c1p[0,-1,1], c2p[0,-1,1], c3p[0,-1,1]})",
            "Chebychev::backgroundFail(mass, {c1f[0,-1,1], c2f[0,-1,1], c3f[0,-1,1]})",
            "efficiency[0.9,0,1]",
            "signalFractionInPassing[0.5]"
        )

    ),

    binnedFit = cms.bool(True),
    binsForFit = cms.uint32(40),
)

matches = [ (1.0,0.4), (0.7,0.3), (0.5,0.2), (0.3,0.15), (0.1,0.1) ]
effs    = [ "", "_NoJPsi", "_NoBestJPsi" ]
if False:
    matches = [ (1.0,0.4), (0.7,0.3), (0.5,0.2), (0.3,0.15), (0.1,0.1) ]
if False:
    effs   =  [ "_NoBestJPsi" ]
    matches = [ (0.3,0.15) ]

process.TnP_Tracking = Template.clone(
    InputFileNames = cms.vstring(),
    InputDirectoryName = cms.string("tpTreeSta"),
    InputTreeName = cms.string("fitter_tree"),
    OutputFileName = cms.string("TnP_Tracking_%s.root" % scenario),
    Efficiencies = cms.PSet()
)
if noEta: process.TnP_Tracking.OutputFileName = "TnP_Tracking_NoEta_%s.root" % scenario

PREFIX="root://pcmssd12.cern.ch//data/gpetrucc/7TeV/tnp/JPsi-2010.01.25/"
PREFIX="/afs/cern.ch/user/g/gpetrucc/scratch0/tnp/CMSSW_3_9_7/src/MuonAnalysis/TagAndProbe/test/jpsi/"
PREFIX="/data/gpetrucc/7TeV/tnp/JPsi-2010.01.25/"
PREFIX="root://castorpublic.cern.ch//castor/cern.ch/user/g/gpetrucc/TnP/JPsiMuMu/"
if scenario == "data_all":
    process.TnP_Tracking.InputFileNames = cms.vstring(
        PREFIX+'tnpJPsi_Data2010B_Nov4.root',
    )
elif scenario == "signal_mc":
    process.TnP_Tracking.InputFileNames = [
        PREFIX+"tnpJPsi_JPsiToMuMu_Fall10.root",
    ]
elif scenario == "relval_mc":
    process.TnP_Tracking.InputFileNames = [
        PREFIX+"tnpJPsi_MC_RelVal.root",
    ]


sampleToPdfMap = { "": "gaussPlusCubic", "NoJPsi":"gaussPlusFloatCubic", "NoBestJPsi":"gaussPlusFloatCubic"}
for (dr,de) in matches:
    for X in effs:
        label = "dr%03de%03d%s" % (100*dr, 100*de, X.replace("_",""))
        if len(args) > 1 and args[1] != label: 
            print "Skipping "+label
            continue
        module = process.TnP_Tracking.clone(OutputFileName = cms.string("TnP_Tracking_%s_%s.root" % (label,scenario)))
        setattr(module.Cuts, "dr_cut",   cms.vstring("tk_deltaR"  +X, "tk_deltaR"  +X, str(dr)))
        setattr(module.Cuts, "deta_cut", cms.vstring("tk_deltaEta"+X, "tk_deltaEta"+X, str(de)))
        common = cms.PSet(
            EfficiencyCategoryAndState = cms.vstring("dr_cut","below", "deta_cut","below"),
            UnbinnedVariables = cms.vstring("mass"),
            BinToPDFmap = cms.vstring(sampleToPdfMap[X.replace("_","")])
        )
        setattr(module.Efficiencies, "eff_"    +label, cms.PSet(common, BinnedVariables = ONE_BIN))
        if ichep:
            setattr(module.Efficiencies, "eff_ichep_"    +label, cms.PSet(common, BinnedVariables = ICHEP_BINS))
        if noEta == False:
            setattr(module.Efficiencies, "eff_abseta_"+label, cms.PSet(common, BinnedVariables = ETA_BINS))
            #setattr(module.Efficiencies, "eff_etaphi_"+label, cms.PSet(common, BinnedVariables = ETA_PHI_BINS))
        if False and scenario == "data_all":
            setattr(module.Efficiencies, "eff_run_"+label, cms.PSet(common, BinnedVariables = RUN_BINS))
        if True and scenario == "data_all":
            setattr(module.Efficiencies, "eff_vtx_"+label, cms.PSet(common, BinnedVariables = VTX_BINS))
        setattr(process,"TnP_Tracking_"+label, module)
        setattr(process,"p_TnP_Tracking_"+label, cms.Path(module))

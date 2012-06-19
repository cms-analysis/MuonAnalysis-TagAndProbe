import FWCore.ParameterSet.Config as cms

import sys
args = sys.argv[1:]
if (sys.argv[0] == "cmsRun"): args =sys.argv[2:]
scenario = "data_all"
if len(args) > 0: scenario = args[0]
print "Will run scenario ", scenario 
doEta = True
doVtx = True

CONSTRAINTS = cms.PSet(
    outerValidHits = cms.vstring("pass"),
    tag_IsoMu24_eta2p1 = cms.vstring("pass"),
    tag_combRelIso = cms.vdouble(-1,0.15),
)
ONE_BIN = cms.PSet(CONSTRAINTS,
    pt = cms.vdouble( 10, 120 ),
    eta = cms.vdouble(-2.4, 2.4),
)
ETA_BINS = cms.PSet(CONSTRAINTS,
    pt = cms.vdouble( 10, 120 ),
    #eta = cms.vdouble(-2.4,-2.1,-1.6,-1.1,-0.6, 0, 0.6, 1.1, 1.6, 2.1, 2.4)
    eta = cms.vdouble(*[-2.4+0.1*x for x in range(0,49)])
)
VTX_BINS = ONE_BIN.clone(
    tag_nVertices = cms.vdouble(*[2*i+0.5 for i in xrange(21)])
)


process = cms.Process("TagProbe")

process.load('FWCore.MessageService.MessageLogger_cfi')

process.source = cms.Source("EmptySource")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

Template = cms.EDAnalyzer("TagProbeFitTreeAnalyzer",
    NumCPU = cms.uint32(1),
    SaveWorkspace = cms.bool(False),

    Variables = cms.PSet(
        mass = cms.vstring("Tag-Probe Mass", "60", "130", "GeV/c^{2}"),
        pt     = cms.vstring("Probe p_{T}", "0", "1000", "GeV/c"),
        p      = cms.vstring("Probe p", "0", "1000", "GeV/c"),
        eta    = cms.vstring("Probe #eta", "-2.5", "2.5", ""),
        abseta = cms.vstring("Probe |#eta|", "0", "2.5", ""),
        phi    = cms.vstring("Probe #phi", "-3.1416", "3.1416", ""),
        tag_pt = cms.vstring("Tag p_{T}", "2.6", "1000", "GeV/c"),
        tk_deltaR   = cms.vstring("Match #Delta R",    "0", "1000", ""),
        tk_deltaEta = cms.vstring("Match #Delta #eta", "0", "1000", ""),
        tk_deltaR_NoZ   = cms.vstring("Unmatch #Delta R",    "0", "1000", ""),
        tk_deltaEta_NoZ = cms.vstring("Unmatch #Delta #eta", "0", "1000", ""),
        tag_nVertices = cms.vstring("Number of vertices", "0", "999", ""),
        tag_combRelIso = cms.vstring("Number of vertices", "-1", "999", ""),
        run    = cms.vstring("Run Number", "132440", "999999", ""),
    ),

    Categories = cms.PSet(
        tag_IsoMu24_eta2p1 = cms.vstring("hasValidHits",  "dummy[pass=1,fail=0]"),
        outerValidHits = cms.vstring("hasValidHits",  "dummy[pass=1,fail=0]"),
        mcTrue = cms.vstring("MC true", "dummy[true=1,false=0]"),
    ),

    Expressions = cms.PSet(),
    Cuts = cms.PSet(),

    PDFs = cms.PSet(
        voigtPlusExpo = cms.vstring(
            "Voigtian::signal(mass, meanP[90,80,100], width[2.495], sigma[5,1,12])",
            "Exponential::backgroundPass(mass, lp[0,-5,5])",
            "Exponential::backgroundFail(mass, lf[0,-5,5])",
            "efficiency[0.9,0,1]",
            "signalFractionInPassing[0.9]"
        ),
    ),

    binnedFit = cms.bool(True),
    binsForFit = cms.uint32(20),
    saveDistributionsPlot = cms.bool(False),
)

matches = [ (0.3, 0.3) ]
effs    = [ "", "_NoZ" ]

process.TnP_Tracking = Template.clone(
    InputFileNames = cms.vstring(),
    InputDirectoryName = cms.string("tpTreeSta"),
    InputTreeName = cms.string("fitter_tree"),
    OutputFileName = cms.string("TnP_Tracking_%s.root" % scenario),
    Efficiencies = cms.PSet()
)

PREFIX="root://pcmssd12//data/gpetrucc/8TeV/tnp/All/"
if "data_all" in scenario:
    process.TnP_Tracking.InputFileNames = cms.vstring(
        PREFIX+'tnpZ_Data_Run2011A.root',
        PREFIX+'tnpZ_Data_Run2011B_upTo194479.root',
        PREFIX+'tnpZ_Data_Run2011B_194480-195016.root',
    )
if "signal_mc" in scenario:
    process.TnP_Tracking.InputFileNames = [
        PREFIX+"tnpZ_MC_DYJetsToLL_madgraph_Summer12_PU_S7_START52_V9_v2.root",
    ]

sampleToPdfMap = { "": "voigtPlusExpo", "NoZ":"voigtPlusExpo"}
for (dr,de) in matches:
    for X in effs:
        label = "dr%03de%03d%s" % (100*dr, 100*de, X.replace("_",""))
        module = process.TnP_Tracking.clone(OutputFileName = cms.string("TnP_Tracking_%s_%s.root" % (label,scenario)))
        setattr(module.Cuts, "dr_cut",   cms.vstring("tk_deltaR"  +X, "tk_deltaR"  +X, str(dr)))
        common = cms.PSet(
            EfficiencyCategoryAndState = cms.vstring("dr_cut","below"),
            UnbinnedVariables = cms.vstring("mass"),
            BinToPDFmap = cms.vstring(sampleToPdfMap[X.replace("_","")])
        )
        if de < dr:
            setattr(module.Cuts, "deta_cut", cms.vstring("tk_deltaEta"+X, "tk_deltaEta"+X, str(de)))
            common.EfficiencyCategoryAndState += [ "deta_cut","below" ]
        if "avg" in scenario: setattr(module.Efficiencies, "eff_"    +label, cms.PSet(common, BinnedVariables = ONE_BIN))
        if "eta" in scenario: setattr(module.Efficiencies, "eff_eta_"+label, cms.PSet(common, BinnedVariables = ETA_BINS))
        if "vtx" in scenario: setattr(module.Efficiencies, "eff_vtx_"+label, cms.PSet(common, BinnedVariables = VTX_BINS))
        setattr(process,"TnP_Tracking_"+label, module)
        setattr(process,"p_TnP_Tracking_"+label, cms.Path(module))

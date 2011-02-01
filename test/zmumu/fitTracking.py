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
)
ONE_BIN = cms.PSet(CONSTRAINTS,
    pt = cms.vdouble( 10, 120 ),
    eta = cms.vdouble(-2.4, 2.4),
)
ETA_BINS = ONE_BIN.clone(
   eta = cms.vdouble(-2.4, -2.0, -1.6, -1.2, -0.8, -0.4, 0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4)
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
        run    = cms.vstring("Run Number", "132440", "999999", ""),
    ),

    Categories = cms.PSet(
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

PREFIX="root://pcmssd12.cern.ch//data/gpetrucc/7TeV/tnp/JPsi-2010.01.25/"
PREFIX="/data/gpetrucc/7TeV/tnp/JPsi-2010.01.25/"
PREFIX="root://castorpublic.cern.ch//castor/cern.ch/user/g/gpetrucc/TnP/JPsiMuMu/"
PREFIX="/afs/cern.ch/user/g/gpetrucc/scratch0/tnp/CMSSW_3_9_7/src/MuonAnalysis/TagAndProbe/test/zmumu/"
PREFIX=""
if scenario == "data_all":
    process.TnP_Tracking.InputFileNames = cms.vstring(
        PREFIX+'tnpZ_Data2010B_Nov4.root',
    )
elif scenario == "signal_mc":
    process.TnP_Tracking.InputFileNames = [
        PREFIX+"tnpZ_DYToMMpowhegZ2_Fall10PU.root",
    ]

sampleToPdfMap = { "": "voigtPlusExpo", "NoZ":"voigtPlusExpo"}
for (dr,de) in matches:
    for X in effs:
        label = "dr%03de%03d%s" % (100*dr, 100*de, X.replace("_",""))
        module = process.TnP_Tracking.clone(OutputFileName = cms.string("TnP_Tracking_%s_%s.root" % (label,scenario)))
        setattr(module.Cuts, "dr_cut",   cms.vstring("tk_deltaR"  +X, "tk_deltaR"  +X, str(dr)))
        setattr(module.Cuts, "deta_cut", cms.vstring("tk_deltaEta"+X, "tk_deltaEta"+X, str(de)))
        common = cms.PSet(
            EfficiencyCategoryAndState = cms.vstring("dr_cut","below", "deta_cut","below"),
            UnbinnedVariables = cms.vstring("mass"),
            BinToPDFmap = cms.vstring(sampleToPdfMap[X.replace("_","")])
        )
        #setattr(module.Efficiencies, "eff_"    +label, cms.PSet(common, BinnedVariables = ONE_BIN))
        if doEta:
            setattr(module.Efficiencies, "eff_abseta_"+label, cms.PSet(common, BinnedVariables = ETA_BINS))
            #setattr(module.Efficiencies, "eff_etaphi_"+label, cms.PSet(common, BinnedVariables = ETA_PHI_BINS))
        if False and scenario == "data_all":
            setattr(module.Efficiencies, "eff_run_"+label, cms.PSet(common, BinnedVariables = RUN_BINS))
        if doVtx:
            setattr(module.Efficiencies, "eff_vtx_"+label, cms.PSet(common, BinnedVariables = VTX_BINS))
        setattr(process,"TnP_Tracking_"+label, module)
        setattr(process,"p_TnP_Tracking_"+label, cms.Path(module))

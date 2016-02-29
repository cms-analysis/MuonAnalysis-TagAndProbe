import FWCore.ParameterSet.Config as cms

import sys
args = sys.argv[1:]
if (sys.argv[0] == "cmsRun"): args =sys.argv[2:]
scenario = "data_all"
if len(args) > 0: scenario = args[0]
print "Will run scenario ", scenario 
doEta = False
doIchep = True
doVtx = False
allMatches = True

CONSTRAINTS = cms.PSet(
    outerValidHits = cms.vstring("pass"),
    MuX_L2Mu0_L2     = cms.vstring("pass"),
    tag_Mu5_L2Mu0_Mu = cms.vstring("pass"),
)
ONE_BIN = cms.PSet(CONSTRAINTS,
    pt = cms.vdouble( 0, 20 ),
    eta = cms.vdouble(-2.4, 2.4),
)
FOUR_BINS = cms.PSet(CONSTRAINTS,
    pt = cms.vdouble( 0, 20 ),
    abseta = cms.vdouble(0, 1.0, 1.5, 2.0, 2.4),
)
ETA_BINS = ONE_BIN.clone(
   eta = cms.vdouble(-2.4, -2.1, -1.8, -1.5, -1.2, -0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9, 1.2, 1.5, 1.8, 2.1, 2.4)
)
#ICHEP_BINS = cms.PSet(CONSTRAINTS,
#    pt = cms.vdouble( 0, 20 ),
#    abseta = cms.vdouble(0, 0.6, 1.1, 1.6, 2.1, 2.4)
#)
AETA_BINS = cms.PSet(CONSTRAINTS,
    pt = cms.vdouble( 0, 20 ),
    #abseta = cms.vdouble(0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4),
    abseta = cms.vdouble(0, 0.6, 1.1, 1.6, 2.1, 2.4),
)

ETA2_BINS = cms.PSet(CONSTRAINTS,
    pt = cms.vdouble( 0, 20 ),
    eta = cms.vdouble(-2.4,-2.1,-1.6,-1.1,-0.6, 0, 0.6, 1.1, 1.6, 2.1, 2.4),
)
ETA3_BINS = cms.PSet(CONSTRAINTS,
    pt = cms.vdouble( 0, 20 ),
    eta = cms.vdouble(-2.4, -2.1, -1.6, -1.2, -0.9, -0.6, -0.3, -0.2, 0.2, 0.3, 0.6, 0.9, 1.2, 1.6, 2.1, 2.4),
)

ETA_PHI_BINS = cms.PSet(CONSTRAINTS,
    pt = cms.vdouble( 0, 20 ),
#    phi = cms.vdouble(*[-3.14+6.28*x/8 for x in range(0,9)]),
    phi = cms.vdouble(-3.14, -2.5, -2, -1.25, -0.63, 0, 0.63, 1.25, 2., 2.5, 3.14),
    eta = cms.vdouble(-2.4,-0.5,0.5,2.4)
)

VTX_BINS = ONE_BIN.clone(
    tag_nVertices = cms.vdouble(*[2*i+0.5 for i in xrange(15)])
#    tag_nVertices = cms.vdouble(0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5)
#    tag_nVertices = cms.vdouble(9.5, 11.5, 13.5, 15.5, 17.5, 19.5, 21.5, 23.5, 25.5, 27.5, 29.5)
)

PHI_BINS = ONE_BIN.clone(
    phi = cms.vdouble(-3.14, -2.5, -2, -1.25, -0.63, 0, 0.63, 1.25, 2., 2.5, 3.14)
)
PT_BINS = cms.PSet(CONSTRAINTS,
    pt = cms.vdouble(  3, 5, 8, 10, 12, 15, 20 ),
    eta = cms.vdouble(-2.4, 2.4),
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
        tk0_deltaR   = cms.vstring("Match #Delta R",    "0", "1000", ""),
        tk0_deltaEta = cms.vstring("Match #Delta #eta", "0", "1000", ""),
        tk0_deltaR_NoJPsi   = cms.vstring("Unmatch #Delta R",    "0", "1000", ""),
        tk0_deltaEta_NoJPsi = cms.vstring("Unmatch #Delta #eta", "0", "1000", ""),
        tk0_deltaR_NoBestJPsi   = cms.vstring("Unmatch #Delta R",    "0", "1000", ""),
        tk0_deltaEta_NoBestJPsi   = cms.vstring("Unmatch #Delta R",    "0", "1000", ""),
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
    saveDistributionsPlot = cms.bool(False),

)

matches = [ (0.3,0.15) ]
effs    = [ "", "_NoJPsi", "_NoBestJPsi" ]
tracks  = [ "tk", "tk0" ]
if allMatches:
    matches = [ (1.0,0.4), (0.7,0.3), (0.5,0.2), (0.3,0.15), (0.1,0.1) ]

process.TnP_Tracking = Template.clone(
    InputFileNames = cms.vstring(),
    InputDirectoryName = cms.string("tpTreeSta"),
    InputTreeName = cms.string("fitter_tree"),
    OutputFileName = cms.string("TnP_Tracking_%s.root" % scenario),
    Efficiencies = cms.PSet()
)

PREFIX="root://lxcms06//data2/b/botta/TnP_Paper2010_Tree/"
PREFIX="/data/gpetrucc/7TeV/tnp/2011.02.17/"
if scenario == "data_all":
    process.TnP_Tracking.InputFileNames = cms.vstring(
        PREFIX+'tnpJPsi_Data_Nov4B.root',
    )
elif scenario == "signal_mc":
    process.TnP_Tracking.InputFileNames = [
        PREFIX+"tnpJPsi_MC_Prompt.crop.root",
    ]
elif scenario == "beauty_mc":
    process.TnP_Tracking.InputFileNames = [
        PREFIX+"tnpJPsi_MC_Bp.root",
    ]
if "TestData_" in scenario:
    process.TnP_Tracking.InputFileNames = [ "tnpJPsi_Data_%s.root" % ("_".join(scenario.replace("TestData_","").split("_")[:-1]),) ]
if "TestMC_" in scenario:
    process.TnP_Tracking.InputFileNames = [ "tnpJPsi_MC_%s.root" % ("_".join(scenario.replace("TestMC_","").split("_")[:-1]),) ]
if "TestMCGen_" in scenario:
    process.TnP_Tracking.InputFileNames = [ "tnpZ_MC_%s.root" % ("_".join(scenario.replace("TestMCGen_","").split("_")[:-1]),) ]
    process.TnP_Tracking.InputDirectoryName = "tpTreeGen"
    del process.TnP_Tracking.Categories.outerValidHits
    del process.TnP_Tracking.Variables.staValidStations
if "TestMCL1_" in scenario:
    process.TnP_Tracking.InputFileNames = [ "tnpZ_MC_%s.root" % ("_".join(scenario.replace("TestMCL1_","").split("_")[:-1]),) ]
    process.TnP_Tracking.InputDirectoryName = "tpTreeL1"
    del process.TnP_Tracking.Categories.outerValidHits
    del process.TnP_Tracking.Variables.staValidStations



sampleToPdfMap = { "": "gaussPlusCubic", "NoJPsi":"gaussPlusFloatCubic", "NoBestJPsi":"gaussPlusFloatCubic"}
for (dr,de) in matches:
  for tk in tracks:
    for X in effs:
        label = "dr%03de%03d%s" % (100*dr, 100*de, X.replace("_",""))
        if tk != "tk": label  = tk+"_"+label
        if len(args) > 1 and args[1] not in label: continue
        module = process.TnP_Tracking.clone(OutputFileName = cms.string("TnP_Tracking_%s_%s.root" % (label,scenario)))
        setattr(module.Cuts, "dr_cut",   cms.vstring(tk+"_deltaR"  +X, tk+"_deltaR"  +X, str(dr)))
        common = cms.PSet(
            EfficiencyCategoryAndState = cms.vstring("dr_cut","below"),
            UnbinnedVariables = cms.vstring("mass"),
            BinToPDFmap = cms.vstring(sampleToPdfMap[X.replace("_","")])
        )
        if "signal_mc" in scenario: common.UnbinnedVariables.append("weight")
        if de < dr:
            setattr(module.Cuts, "deta_cut", cms.vstring(tk+"_deltaEta"+X, tk+"_deltaEta"+X, str(de)))
            common.EfficiencyCategoryAndState += [ "deta_cut","below" ]
        #setattr(module.Efficiencies, "eff_"    +label, cms.PSet(common, BinnedVariables = ONE_BIN))
        if "ichep" in scenario:
            setattr(module.Efficiencies, "eff_ichep_"    +label, cms.PSet(common, BinnedVariables = ICHEP_BINS))
        if "et1" in scenario:
            setattr(module.Efficiencies, "eff_eta_"+label, cms.PSet(common, BinnedVariables = ETA_BINS))
        if "et2" in scenario:
            setattr(module.Efficiencies, "eff_eta2_"+label, cms.PSet(common, BinnedVariables = ETA2_BINS))
        if "et3" in scenario:
            setattr(module.Efficiencies, "eff_eta3_"+label, cms.PSet(common, BinnedVariables = ETA3_BINS))
        if "aeta" in scenario:
            setattr(module.Efficiencies, "eff_abseta_"+label, cms.PSet(common, BinnedVariables = AETA_BINS))
        if "vtx" in scenario:
            setattr(module.Efficiencies, "eff_vtx_"+label, cms.PSet(common, BinnedVariables = VTX_BINS))
        if "phi" in scenario: 
            setattr(module.Efficiencies, "eff_phi_"+label, cms.PSet(common, BinnedVariables = PHI_BINS))
        if "eph" in scenario:
            setattr(module.Efficiencies, "eff_eph_"+label, cms.PSet(common, BinnedVariables = ETA_PHI_BINS))
        if "pt" in scenario:
            setattr(module.Efficiencies, "eff_pt_"+label, cms.PSet(common, BinnedVariables = PT_BINS))
        if "four" in scenario: 
	  setattr(module.Efficiencies, "eff_four_"+label, cms.PSet(common, BinnedVariables = FOUR_BINS))
        setattr(process,"TnP_Tracking_"+label, module)
        setattr(process,"p_TnP_Tracking_"+label, cms.Path(module))
        if module.InputDirectoryName.value() == "tpTreeGen":
            for l in module.Efficiencies.parameterNames_():
                effpset = getattr(module.Efficiencies,l)
                if not hasattr(getattr(module.Efficiencies,l), 'BinnedVariables'): continue
                if hasattr(getattr(module.Efficiencies,l).BinnedVariables, 'outerValidHits'):
                    del getattr(module.Efficiencies,l).BinnedVariables.outerValidHits
                if hasattr(getattr(module.Efficiencies,l).BinnedVariables, 'staValidStations'):
                    del getattr(module.Efficiencies,l).BinnedVariables.staValidStations
        if module.InputDirectoryName.value() == "tpTreeL1":
            for l in module.Efficiencies.parameterNames_():
                effpset = getattr(module.Efficiencies,l)
                if not hasattr(getattr(module.Efficiencies,l), 'BinnedVariables'): continue
                if hasattr(getattr(module.Efficiencies,l).BinnedVariables, 'outerValidHits'):
                    del getattr(module.Efficiencies,l).BinnedVariables.outerValidHits
                if hasattr(getattr(module.Efficiencies,l).BinnedVariables, 'staValidStations'):
                    del getattr(module.Efficiencies,l).BinnedVariables.staValidStations


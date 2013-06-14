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
    #tag_IsoMu24_eta2p1 = cms.vstring("pass"),
    #tag_combRelIso = cms.vdouble(-1,0.15),
)
ONE_BIN = cms.PSet(CONSTRAINTS,
    pt = cms.vdouble( 10, 120 ),
    eta = cms.vdouble(-2.4, 2.4),
)
TWO_BINS = cms.PSet(CONSTRAINTS,
    pt = cms.vdouble( 10, 120 ),
    abseta = cms.vdouble(0, 1.2, 2.4),
)

ETA_BINS = cms.PSet(CONSTRAINTS,
    pt = cms.vdouble( 10, 120 ),
    #eta = cms.vdouble(-2.4,-2.1,-1.6,-1.1,-0.6, 0, 0.6, 1.1, 1.6, 2.1, 2.4)
    #eta = cms.vdouble(*[-2.4+0.1*x for x in range(0,49)])
    eta = cms.vdouble(*[-2.4+0.2*x for x in range(0,25)])
)
ETA2_BINS = cms.PSet(CONSTRAINTS,
    pt = cms.vdouble( 10, 120 ),
    eta = cms.vdouble(-2.4,-2.1,-1.6,-1.1,-0.6, 0, 0.6, 1.1, 1.6, 2.1, 2.4),
)
ETA3_BINS = cms.PSet(CONSTRAINTS,
    pt = cms.vdouble( 10, 120 ),
    eta = cms.vdouble(-2.4, -2.1, -1.6, -1.2, -0.9, -0.6, -0.3, -0.2, 0.2, 0.3, 0.6, 0.9, 1.2, 1.6, 2.1, 2.4),
)

ETA_PHI_BINS = cms.PSet(CONSTRAINTS,
    pt = cms.vdouble( 20, 120 ),
    eta = cms.vdouble(-2.4,-1.9,1.9,2.4),
    phi = cms.vdouble(*[-3.14+6.28*x/8 for x in range(0,9)]),
)

VTX_BINS = ONE_BIN.clone(
    tag_nVertices = cms.vdouble(*[2*i+0.5 for i in xrange(15)])
)
VTX_BINS.tag_nVertices += [30.5, 35.5]

RUN_BINS_ALL = cms.PSet(CONSTRAINTS,
    pt = cms.vdouble( 10, 120 ),
    run = cms.vdouble(190000, 191280, 192000, 193200, 194000, 194500, 195500, 196100, 197000),
)
RUN_BINS_ETA_PHI = cms.PSet(CONSTRAINTS,
    pt = cms.vdouble( 20, 120 ),
    eta = cms.vdouble(-2.4,-1.9,1.9,2.4),
    phi = cms.vdouble(-3.14, -1.57, 0, 1.57, 3.14),
    run = cms.vdouble(190000, 192000, 193200, 194500, 195500, 197000),
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
        #tag_pt = cms.vstring("Tag p_{T}", "2.6", "1000", "GeV/c"),
        tk_deltaR   = cms.vstring("Match #Delta R",    "0", "1000", ""),
        tk_deltaEta = cms.vstring("Match #Delta #eta", "0", "1000", ""),
        tk_deltaR_NoZ   = cms.vstring("Unmatch #Delta R",    "0", "1000", ""),
        tk_deltaEta_NoZ = cms.vstring("Unmatch #Delta #eta", "0", "1000", ""),
        tag_nVertices = cms.vstring("Number of vertices", "0", "999", ""),
        #tag_combRelIso = cms.vstring("Number of vertices", "-1", "999", ""),
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

if "two" in scenario:
    matches = [ (0.1,0.05), (0.2,0.1), (0.2,0.2), (0.3, 0.3), (0.6,0.4), (1.0,0.4) ]
    

process.TnP_Tracking = Template.clone(
    InputFileNames = cms.vstring(),
    InputDirectoryName = cms.string("tpTreeSta"),
    InputTreeName = cms.string("fitter_tree"),
    OutputFileName = cms.string("TnP_Tracking_%s.root" % scenario),
    Efficiencies = cms.PSet()
)

PREFIX="root://pcmssd12//data/gpetrucc/8TeV/tnp/Tracking/"
if "data_ABCD" in scenario:
    process.TnP_Tracking.InputFileNames = [ 
        'root://pcmssd12.cern.ch//data/gpetrucc/8TeV/tnp/Tracking/tnpZ_run2012A_13Jul2012_forTracking.root',
        'root://pcmssd12.cern.ch//data/gpetrucc/8TeV/tnp/Tracking/tnpZ_run2012B_13Jul2012_forTracking.root',
        'root://pcmssd12.cern.ch//data/gpetrucc/8TeV/tnp/Tracking/tnpZ_run2012Cv1_24Aug2012_forTracking.root',
        'root://pcmssd12.cern.ch//data/gpetrucc/8TeV/tnp/Tracking/tnpZ_run2012Cv2_upTo203002_forTracking.root',
        'root://pcmssd12.cern.ch//data/gpetrucc/8TeV/tnp/Tracking/tnpZ_run2012D_upTo207898_forTracking.root',
        'root://pcmssd12.cern.ch//data/gpetrucc/8TeV/tnp/Tracking/tnpZ_run2012D_207899-208686_forTracking.root',
    ]
if "data_2012C" in scenario:
    process.TnP_Tracking.InputFileNames = [ PREFIX+'tnpZ_Data_Run2012C_forTracking.root', ]
if "data_2012D" in scenario:
    process.TnP_Tracking.InputFileNames = [ 
        'root://pcmssd12.cern.ch//data/gpetrucc/8TeV/tnp/Tracking/tnpZ_run2012D_upTo207898_forTracking.root',
        'root://pcmssd12.cern.ch//data/gpetrucc/8TeV/tnp/Tracking/tnpZ_run2012D_207899-208686_forTracking.root',
    ]
if "signal_mc" in scenario:
    process.TnP_Tracking.InputFileNames = [
        'root://pcmssd12.cern.ch//data/gpetrucc/8TeV/tnp/Tracking/tnpZ_MC53X_DYJets_forTracking_withNVtxWeightsUpTo203002.root' 
    ]
    process.TnP_Tracking.WeightVariable = cms.string("weight")
    process.TnP_Tracking.Variables.weight = cms.vstring("weight","0","10","")
if "signal_mc_ABCD" in scenario:
    process.TnP_Tracking.InputFileNames = [
        'root://pcmssd12.cern.ch//data/gpetrucc/8TeV/tnp/Tracking/tnpZ_MC53X_DYJets_powheg_forTracking_withNVtxWeightsABCD.root' 
    ]
    process.TnP_Tracking.WeightVariable = cms.string("weight")
    process.TnP_Tracking.Variables.weight = cms.vstring("weight","0","10","")
if "signal_mc_2012D" in scenario:
    process.TnP_Tracking.InputFileNames = [
        'root://pcmssd12.cern.ch//data/gpetrucc/8TeV/tnp/Tracking/tnpZ_MC53X_DYJets_forTracking_withNVtxWeights2012D.root' 
    ]
    process.TnP_Tracking.WeightVariable = cms.string("weight")
    process.TnP_Tracking.Variables.weight = cms.vstring("weight","0","10","")


if "muReRec2012Cv1R" in scenario:
    process.TnP_Tracking.InputFileNames = [ "root://eoscms//eos/cms/store/caf/user/gpetrucc/tnpZ_Data_R.root" ]
if "muReRec2012Cv1P" in scenario:
    process.TnP_Tracking.InputFileNames = [ "root://eoscms//eos/cms/store/caf/user/gpetrucc/tnpZ_Data_P.root" ]



sampleToPdfMap = { "": "voigtPlusExpo", "NoZ":"voigtPlusExpo"}
for (dr,de) in matches:
    for X in effs:
        label = "dr%03de%03d%s" % (100*dr, 100*de, X.replace("_",""))
        if len(args) > 1 and args[1] not in label: continue
        module = process.TnP_Tracking.clone(OutputFileName = cms.string("TnP_Tracking_%s_%s.root" % (label,scenario)))
        setattr(module.Cuts, "dr_cut",   cms.vstring("tk_deltaR"  +X, "tk_deltaR"  +X, str(dr)))
        common = cms.PSet(
            EfficiencyCategoryAndState = cms.vstring("dr_cut","below"),
            UnbinnedVariables = cms.vstring("mass"),
            BinToPDFmap = cms.vstring(sampleToPdfMap[X.replace("_","")])
        )
        if "signal_mc" in scenario: common.UnbinnedVariables.append("weight")
        if de < dr:
            setattr(module.Cuts, "deta_cut", cms.vstring("tk_deltaEta"+X, "tk_deltaEta"+X, str(de)))
            common.EfficiencyCategoryAndState += [ "deta_cut","below" ]
        if "avg" in scenario: setattr(module.Efficiencies, "eff_"    +label, cms.PSet(common, BinnedVariables = ONE_BIN))
        if "eta" in scenario: setattr(module.Efficiencies, "eff_eta_"+label, cms.PSet(common, BinnedVariables = ETA_BINS))
        if "et2" in scenario: setattr(module.Efficiencies, "eff_eta2_"+label, cms.PSet(common, BinnedVariables = ETA2_BINS))
        if "et3" in scenario: setattr(module.Efficiencies, "eff_eta3_"+label, cms.PSet(common, BinnedVariables = ETA3_BINS))
        if "vtx" in scenario: setattr(module.Efficiencies, "eff_vtx_"+label, cms.PSet(common, BinnedVariables = VTX_BINS))
        if "two" in scenario: setattr(module.Efficiencies, "eff_two_"+label, cms.PSet(common, BinnedVariables = TWO_BINS))
        if "eph" in scenario: setattr(module.Efficiencies, "eff_eph_"+label, cms.PSet(common, BinnedVariables = ETA_PHI_BINS))
        if "run" in scenario: setattr(module.Efficiencies, "eff_run_"+label, cms.PSet(common, BinnedVariables = RUN_BINS_ALL))
        if "rep" in scenario: setattr(module.Efficiencies, "eff_rep_"+label, cms.PSet(common, BinnedVariables = RUN_BINS_ETA_PHI))
        setattr(process,"TnP_Tracking_"+label, module)
        setattr(process,"p_TnP_Tracking_"+label, cms.Path(module))

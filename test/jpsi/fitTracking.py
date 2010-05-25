import FWCore.ParameterSet.Config as cms

import sys
args = sys.argv[1:]
if (sys.argv[0] == "cmsRun"): args =sys.argv[2:]
scenario = "signal_0.1pb"
if len(args) > 0: scenario = args[0]
print "Will run scenario ", scenario 


CONSTRAINTS = cms.PSet(
    hasValidHits = cms.vstring("pass"),
    tag_HLTMu3   = cms.vstring("pass"),
)
PT_ETA_BINS = cms.PSet(
    CONSTRAINTS,
    pt = cms.vdouble( 2, 4.5, 15),
    eta = cms.vdouble(-2.5, -1.1, 1.1, 2.5)
)
ETA_PHI_BINS = cms.PSet(
    CONSTRAINTS,
    phi = cms.vdouble(*[3.1416*i/3.0 for i in range(-3,4)]), # 6 bins
    eta = cms.vdouble(-2.5, -1.1, 1.1, 2.5)
)



process = cms.Process("TagProbe")

process.load('FWCore.MessageService.MessageLogger_cfi')

process.source = cms.Source("EmptySource")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

Template = cms.EDAnalyzer("TagProbeFitTreeAnalyzer",
    NumCPU = cms.uint32(1),
    SaveWorkspace = cms.bool(False),

    Variables = cms.PSet(
        mass = cms.vstring("Tag-Probe Mass", "2.0", "4.3", "GeV/c^{2}"),
        pt = cms.vstring("Probe p_{T}", "0", "1000", "GeV/c"),
        p = cms.vstring("Probe p", "0", "1000", "GeV/c"),
        eta = cms.vstring("Probe #eta", "-2.5", "2.5", ""),
        phi = cms.vstring("Probe #phi", "-3.1416", "3.1416", ""),
        tag_pt = cms.vstring("Tag p_{T}", "2.6", "1000", "GeV/c"),
    ),

    Categories = cms.PSet(
        mcTrue = cms.vstring("MC true", "dummy[true=1,false=0]"),
        passing = cms.vstring("passing", "dummy[pass=1,fail=0]"),
        hasValidHits = cms.vstring("hasValidHits",  "dummy[pass=1,fail=0]"),
        tag_HLTMu3 = cms.vstring("tag_HLTMu3", "dummy[pass=1,fail=0]"),
    ),

    PDFs = cms.PSet(
        gaussPlusCubic = cms.vstring(
            "Gaussian::signal(mass, mean[3.1,3.0,3.2], sigma[0.2,0.05,0.5])",
            "Chebychev::backgroundPass(mass, {c1[0,-1,1], c2[0,-1,1], c3[0,-1,1]})",
            "Chebychev::backgroundFail(mass, {c1,c2,c3})",
            "efficiency[0.9,0,1]",
            "signalFractionInPassing[0.9]"
        )
    )
)

process.TnP_Tracking = Template.clone(
    InputFileNames = cms.vstring("tnpJPsi_JPsiMuMu_Spring10_0.1pb.root"),
    InputDirectoryName = cms.string("histoTracking"),
    InputTreeName = cms.string("fitter_tree"),
    OutputFileName = cms.string("TnP_Tracking_%s.root" % scenario),
    Efficiencies = cms.PSet(
        pt_eta = cms.PSet(
            EfficiencyCategoryAndState = cms.vstring("passing","pass"),
            UnbinnedVariables = cms.vstring("mass"),
            BinnedVariables = PT_ETA_BINS,
            BinToPDFmap = cms.vstring("gaussPlusCubic")
        ),
        pt_eta_mcTrue = cms.PSet(
            EfficiencyCategoryAndState = cms.vstring("passing","pass"),
            UnbinnedVariables = cms.vstring("mass"),
            BinnedVariables = PT_ETA_BINS.clone(
                mcTrue = cms.vstring("true")
            )
        ),
        eta_phi = cms.PSet(
            EfficiencyCategoryAndState = cms.vstring("passing","pass"),
            UnbinnedVariables = cms.vstring("mass"),
            BinnedVariables = ETA_PHI_BINS,
            BinToPDFmap = cms.vstring("gaussPlusCubic")
        ),
        eta_phi_mcTrue = cms.PSet(
            EfficiencyCategoryAndState = cms.vstring("passing","pass"),
            UnbinnedVariables = cms.vstring("mass"),
            BinnedVariables = ETA_PHI_BINS.clone(
                mcTrue = cms.vstring("true")
            )
        ),
    )
)
if scenario == "all_0.1pb":
    process.TnP_Tracking.InputFileNames = [ 
        "tnpJPsi_JPsiMuMu_Spring10_0.1pb.root",
        "tnpJPsi_ppMuX_Spring10_0.1pb.root"
    ]
elif scenario == "signal_0.5pb":
    process.TnP_Tracking.InputFileNames = [ "tnpJPsi_JPsiMuMu_Spring10_0.5pb.root" ]

if scenario.startswith("data"):
    if scenario == "data_all":
        process.TnP_Tracking.InputFileNames = [ "tnpJPsi_Data.root" ]
    elif scenario == "data_0.001pb":
         process.TnP_Tracking.InputFileNames = [ "/afs/cern.ch/user/g/gpetrucc/scratch0/tnp/tnpJPsi_Data_fromMay6th.root" ]
    elif scenario == "datalike_mc":
         process.TnP_Tracking.InputFileNames = [ "tnpJPsi_JPsiMuMu_Spring10_0.03pb.root", "tnpJPsi_ppMuX_Spring10_0.03pb.root" ]
    process.TnP_Tracking.Variables.tag_pt[1] = "0.0"; # don't cut on tag pt
    process.TnP_Tracking.Efficiencies = cms.PSet(
        eff = cms.PSet(
            EfficiencyCategoryAndState = cms.vstring("passing","pass"),
            UnbinnedVariables = cms.vstring("mass"),
            BinnedVariables = cms.PSet(
                pt  = cms.vdouble(3,15),
                eta = cms.vdouble(-2.5,2.5),
                hasValidHits = cms.vstring("pass"),
            ),
            BinToPDFmap = cms.vstring("gaussPlusCubic")
        ),
        eff_eta = cms.PSet(
            EfficiencyCategoryAndState = cms.vstring("passing","pass"),
            UnbinnedVariables = cms.vstring("mass"),
            BinnedVariables = cms.PSet(
                pt  = cms.vdouble(3,15),
                eta = cms.vdouble(-2.5,-1.1,1.1,2.5),
                hasValidHits = cms.vstring("pass"),
            ),
            BinToPDFmap = cms.vstring("gaussPlusCubic")
        )
    )

process.p = cms.Path(
    process.TnP_Tracking
)


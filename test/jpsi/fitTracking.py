import FWCore.ParameterSet.Config as cms

FILEPREFIX = "signal_"
CONSTRAINTS = cms.PSet(
    hasValidHits = cms.vstring("pass"),
    tag_HLTMu3   = cms.vstring("pass"),
)
PT_ETA_BINS = cms.PSet(
    CONSTRAINTS,
    pt = cms.vdouble( 2, 4.5, 20),
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
        #tag_L1DiMuOpen = cms.vstring("tag_L1DoubleMuOpen", "dummy[pass=1,fail=0]")
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
    InputFileNames = cms.vstring(
        "tnpJPsi_JPsiMuMu_Spring10_0.5pb.root",
        #"tnpJPsi_ppMuX_Spring10_0.1pb.root"
    ),
    InputDirectoryName = cms.string("histoTracking"),
    InputTreeName = cms.string("fitter_tree"),
    OutputFileName = cms.string(FILEPREFIX+"TnP_Tracking_0.5pb.root"),
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

process.p = cms.Path(
    process.TnP_Tracking
)


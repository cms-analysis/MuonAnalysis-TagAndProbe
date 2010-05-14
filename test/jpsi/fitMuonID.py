import FWCore.ParameterSet.Config as cms

FILEPREFIX = "signal_"
#FILEPREFIX = "withbg_"
CONSTRAINTS = cms.PSet(
    #Glb        = cms.vstring("true"),
    tag_HLTMu3 = cms.vstring("pass"),
)
PT_ETA_BINS = cms.PSet(
    CONSTRAINTS,
    pt = cms.vdouble( 2, 3, 4.5, 6, 20),
    eta = cms.vdouble(-2.4, -1.1, 1.1, 2.4)
)

process = cms.Process("TagProbe")

process.load('FWCore.MessageService.MessageLogger_cfi')

process.source = cms.Source("EmptySource")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

Template = cms.EDAnalyzer("TagProbeFitTreeAnalyzer",
    NumCPU = cms.uint32(1),
    SaveWorkspace = cms.bool(False),

    Variables = cms.PSet(
        mass = cms.vstring("Tag-Probe Mass", "2.8", "3.5", "GeV/c^{2}"),
        pt = cms.vstring("Probe p_{T}", "0", "1000", "GeV/c"),
        p = cms.vstring("Probe p", "0", "1000", "GeV/c"),
        eta = cms.vstring("Probe #eta", "-2.5", "2.5", ""),
        phi = cms.vstring("Probe #phi", "-3.1416", "3.1416", ""),
        tag_pt = cms.vstring("Tag p_{T}", "2.6", "1000", "GeV/c"),
    ),

    Categories = cms.PSet(
        mcTrue = cms.vstring("MC true", "dummy[true=1,false=0]"),
        POG_TMA = cms.vstring("POG_TMA", "dummy[pass=1,fail=0]"),
        POG_Glb = cms.vstring("POG_Glb", "dummy[pass=1,fail=0]"),
        POG_Sta = cms.vstring("POG_Sta", "dummy[pass=1,fail=0]"),
        tag_HLTMu3 = cms.vstring("tag_HLTMu3", "dummy[pass=1,fail=0]"),
    ),

    PDFs = cms.PSet(
        gaussPlusExpo = cms.vstring(
            "Gaussian::signal(mass, mean[3.1,3.0,3.2], sigma[0.05,0.02,0.1])",
            "Exponential::backgroundPass(mass, lp[0,-1,1])",
            "Exponential::backgroundFail(mass, lf[0,-1,1])",
            "efficiency[0.9,0,1]",
            "signalFractionInPassing[0.9]"
        )
    )
)

process.TnP_MuonID = Template.clone(
    InputFileNames = cms.vstring(
        "tnpJPsi_JPsiMuMu_Spring10_0.5pb.root",
        #"tnpJPsi_JPsiMuMu_Spring10_0.1pb.root",
        #"tnpJPsi_ppMuX_Spring10_0.1pb.root"
    ),
    InputDirectoryName = cms.string("histoMuFromTk"),
    InputTreeName = cms.string("fitter_tree"),
    OutputFileName = cms.string(FILEPREFIX+"TnP_MuonID_0.5pb.root"),
    #OutputFileName = cms.string(FILEPREFIX+"TnP_MuonID_0.1pb.root"),
    Efficiencies = cms.PSet(),
)

for T in [ "POG_Glb", "POG_TMA", "POG_Sta" ]:
    setattr(process.TnP_MuonID.Efficiencies, T+"_pt_eta", cms.PSet(
            EfficiencyCategoryAndState = cms.vstring(T,"pass"),
            UnbinnedVariables = cms.vstring("mass"),
            BinnedVariables = PT_ETA_BINS,
            BinToPDFmap = cms.vstring("gaussPlusExpo")
    ))
    setattr(process.TnP_MuonID.Efficiencies, T+"_pt_eta_mcTrue", cms.PSet(
            EfficiencyCategoryAndState = cms.vstring(T,"pass"),
            UnbinnedVariables = cms.vstring("mass"),
            BinnedVariables = PT_ETA_BINS.clone(
                mcTrue = cms.vstring("true")
            )
    ))
    
if True:
    process.TnP_MuonID.InputFileNames = [ "tnpJPsi_Data.root" ]
    process.TnP_MuonID.OutputFileName = "data_TnP_MuonID_1nb.root"
    process.TnP_MuonID.Efficiencies = cms.PSet()
    process.TnP_MuonID.Variables.tag_pt[1] = "0.0"; # cut on tag pt
    for T in [ "POG_Glb", "POG_TMA", "POG_Sta" ]:
        setattr(process.TnP_MuonID.Efficiencies, T+"_pt_eta", cms.PSet(
            EfficiencyCategoryAndState = cms.vstring(T,"pass"),
            UnbinnedVariables = cms.vstring("mass"),
            BinnedVariables = cms.PSet(
                pt  = cms.vdouble( 2, 3.5, 15 ),
                eta = cms.vdouble(-2.4, -1.1, 1.1, 2.4),
            ),
            BinToPDFmap = cms.vstring("gaussPlusExpo")
        ))
        setattr(process.TnP_MuonID.Efficiencies, T+"_pt", cms.PSet(
            EfficiencyCategoryAndState = cms.vstring(T,"pass"),
            UnbinnedVariables = cms.vstring("mass"),
            BinnedVariables = cms.PSet(
                pt  = cms.vdouble( 2, 3.5, 15 ),
                eta = cms.vdouble(-2.4, 2.4),
            ),
            BinToPDFmap = cms.vstring("gaussPlusExpo")
        ))


process.p = cms.Path(
    process.TnP_MuonID
)


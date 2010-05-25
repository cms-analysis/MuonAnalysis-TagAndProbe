import FWCore.ParameterSet.Config as cms

import sys
args = sys.argv[1:]
if (sys.argv[0] == "cmsRun"): args =sys.argv[2:]
scenario = "signal_0.1pb"
if len(args) > 0: scenario = args[0]
print "Will run scenario ", scenario 

CONSTRAINTS = cms.PSet(
    tag_HLTMu3 = cms.vstring("pass"),
)
PT_ETA_BINS = cms.PSet(
    CONSTRAINTS,
    pt = cms.vdouble( 2, 3, 4.5, 6, 12),
    eta = cms.vdouble(-2.4, -1.1, 1.1, 2.4)
)
if scenario == "signal_0.5pb":
    PT_ETA_BINS.pt = cms.vdouble( 2, 3, 4.5, 6, 10, 20)


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
    InputFileNames = cms.vstring("tnpJPsi_JPsiMuMu_Spring10_0.1pb.root"),
    InputDirectoryName = cms.string("histoMuFromTk"),
    InputTreeName = cms.string("fitter_tree"),
    OutputFileName = cms.string("TnP_MuonID_%s.root" % scenario),
    Efficiencies = cms.PSet(),
)
if scenario == "all_0.1pb":
    process.TnP_MuonID.InputFileNames = [ 
        "tnpJPsi_JPsiMuMu_Spring10_0.1pb.root",
        "tnpJPsi_ppMuX_Spring10_0.1pb.root"
    ]
elif scenario == "signal_0.5pb":
    process.TnP_MuonID.InputFileNames = [ "tnpJPsi_JPsiMuMu_Spring10_0.5pb.root" ]
  


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


if scenario.startswith("data"):
    if scenario == "data_all":
        process.TnP_MuonID.InputFileNames = [ "tnpJPsi_Data.root" ]
    elif scenario == "data_0.001pb":
        process.TnP_MuonID.InputFileNames = [ "/afs/cern.ch/user/g/gpetrucc/scratch0/tnp/tnpJPsi_Data_fromMay6th.root" ]
        process.TnP_MuonID.Variables.tag_pt[1] = "0.0"; # cut on tag pt
    elif scenario == "datalike_mc":
        process.TnP_MuonID.InputFileNames = [
            "tnpJPsi_JPsiMuMu_Spring10_0.1pb.root",
            "tnpJPsi_ppMuX_Spring10_0.1pb.root"
        ]
    process.TnP_MuonID.Efficiencies = cms.PSet()
    for T in [ "POG_Glb", "POG_TMA", "POG_Sta" ]:
        setattr(process.TnP_MuonID.Efficiencies, T+"_pt_eta", cms.PSet(
            EfficiencyCategoryAndState = cms.vstring(T,"pass"),
            UnbinnedVariables = cms.vstring("mass"),
            BinnedVariables = cms.PSet(
                pt  = cms.vdouble( 2, 3, 5, 12 ),
                eta = cms.vdouble(-2.4, -1.0, 1.0, 2.4),
            ),
            BinToPDFmap = cms.vstring("gaussPlusExpo")
        ))
        setattr(process.TnP_MuonID.Efficiencies, T+"_pt", cms.PSet(
            EfficiencyCategoryAndState = cms.vstring(T,"pass"),
            UnbinnedVariables = cms.vstring("mass"),
            BinnedVariables = cms.PSet(
                pt  = cms.vdouble( 2, 3, 5, 12 ),
                eta = cms.vdouble(-2.4, 2.4),
            ),
            BinToPDFmap = cms.vstring("gaussPlusExpo")
        ))


process.TnP_MuonID_fromCalo = process.TnP_MuonID.clone(
    InputDirectoryName = "histoMuFromCal",
    OutputFileName = "TnP_MuonID_fromCalo_%s.root" % scenario
)
process.p = cms.Path(
    process.TnP_MuonID +
    process.TnP_MuonID_fromCalo
)


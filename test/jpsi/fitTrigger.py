import FWCore.ParameterSet.Config as cms

import sys
args = sys.argv[1:]
if (sys.argv[0] == "cmsRun"): args =sys.argv[2:]
scenario = "signal_0.1pb"
if len(args) > 0: scenario = args[0]
print "Will run scenario ", scenario 

CONSTRAINTS = cms.PSet(
    Glb        = cms.vstring("true"),
    tag_HLTMu3 = cms.vstring("pass"),
)
PT_ETA_BINS = cms.PSet(
    CONSTRAINTS,
    pt = cms.vdouble( 2, 3, 4.5, 6, 12),
    eta = cms.vdouble(-2.1, -1.1, 1.1, 2.1)
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
        Glb    = cms.vstring("GlobalMu","dummy[true=1,false=0]"),
        HLTMu3     = cms.vstring("HLTMu3", "dummy[pass=1,fail=0]"),
        L1DiMuOpen = cms.vstring("HLTMu3", "dummy[pass=1,fail=0]"),
        tag_HLTMu3 = cms.vstring("tag_HLTMu3", "dummy[pass=1,fail=0]"),
    ),

    PDFs = cms.PSet(
        gaussPlusExpo = cms.vstring(
            "Gaussian::signal(mass, mean[3.1,3.0,3.2], sigma[0.05,0.02,0.1])",
            "Exponential::backgroundPass(mass, lp[0,-1,1])",
            "Exponential::backgroundFail(mass, lp)",  # same slope, they're both muons
            #"Exponential::backgroundFail(mass, lf[0,-1,1])",
            "efficiency[0.9,0,1]",
            "signalFractionInPassing[0.9]"
        )
    )
)

process.TnP_Trigger = Template.clone(
    InputFileNames = cms.vstring("tnpJPsi_JPsiMuMu_Spring10_0.1pb.root"),
    InputDirectoryName = cms.string("histoTrigger"),
    InputTreeName = cms.string("fitter_tree"),
    OutputFileName = cms.string("TnP_Trigger_%s.root" % scenario),
    Efficiencies = cms.PSet(),
)
if scenario == "all_0.1pb":
    process.TnP_Trigger.InputFileNames = [ 
        "tnpJPsi_JPsiMuMu_Spring10_0.1pb.root",
        "tnpJPsi_ppMuX_Spring10_0.1pb.root"
    ]
elif scenario == "signal_0.5pb":
    process.TnP_Trigger.InputFileNames = [ "tnpJPsi_JPsiMuMu_Spring10_0.5pb.root" ]
  

for T in [ "HLTMu3", "L1DiMuOpen" ]:
    setattr(process.TnP_Trigger.Efficiencies, T+"_pt_eta", cms.PSet(
            EfficiencyCategoryAndState = cms.vstring(T,"pass"),
            UnbinnedVariables = cms.vstring("mass"),
            BinnedVariables = PT_ETA_BINS,
            BinToPDFmap = cms.vstring("gaussPlusExpo")
    ))
    setattr(process.TnP_Trigger.Efficiencies, T+"_pt_eta_mcTrue", cms.PSet(
            EfficiencyCategoryAndState = cms.vstring(T,"pass"),
            UnbinnedVariables = cms.vstring("mass"),
            BinnedVariables = PT_ETA_BINS.clone(
                mcTrue = cms.vstring("true")
            )
    ))
    
if scenario.startswith("data"):
    if scenario == "data_all":
        process.TnP_Trigger.InputFileNames = [ "tnpJPsi_Data.root" ]
    elif scenario == "data_0.001pb":
         process.TnP_Trigger.InputFileNames = [ "/afs/cern.ch/user/g/gpetrucc/scratch0/tnp/tnpJPsi_Data_fromMay6th.root" ]
    process.TnP_Trigger.Efficiencies = cms.PSet()
    process.TnP_Trigger.Variables.tag_pt[1] = "0.0"; # don't cut on tag pt
    process.TnP_Trigger.Variables.run       = cms.vstring("Run number", "0", "9999999", "");
    for T in [ "HLTMu3", "L1DiMuOpen" ]:
        setattr(process.TnP_Trigger.Efficiencies, T+"_pt_eta", cms.PSet(
            EfficiencyCategoryAndState = cms.vstring(T,"pass"),
            UnbinnedVariables = cms.vstring("mass"),
            BinnedVariables = cms.PSet(
                Glb = cms.vstring("true"),
                pt  = cms.vdouble( 2, 3, 12 ),
                eta = cms.vdouble(-2.1,-1.0,1.0, 2.1),
                run = cms.vdouble(0,9999999),
            ),
            BinToPDFmap = cms.vstring("gaussPlusExpo")
        ))
        setattr(process.TnP_Trigger.Efficiencies, T+"_pt", cms.PSet(
            EfficiencyCategoryAndState = cms.vstring(T,"pass"),
            UnbinnedVariables = cms.vstring("mass"),
            BinnedVariables = cms.PSet(
                Glb = cms.vstring("true"),
                pt  = cms.vdouble( 2, 3, 15 ),
                eta = cms.vdouble(-2.1, 2.1),
                run = cms.vdouble(0,9999999),
            ),
            BinToPDFmap = cms.vstring("gaussPlusExpo")
        ))
        if T.find("HLT") != -1:
            getattr(process.TnP_Trigger.Efficiencies, T+"_pt_eta").BinnedVariables.run = cms.vdouble(133443,9999999)
            getattr(process.TnP_Trigger.Efficiencies, T+"_pt"    ).BinnedVariables.run = cms.vdouble(133443,9999999)
        setattr(process.TnP_Trigger.Efficiencies, T+"_run", cms.PSet(
            EfficiencyCategoryAndState = cms.vstring(T,"pass"),
            UnbinnedVariables = cms.vstring("mass"),
            BinnedVariables = cms.PSet(
                Glb = cms.vstring("true"),
                pt  = cms.vdouble( 3, 15),
                eta = cms.vdouble(-2.1, 2.1),
                run = cms.vdouble(132440,133443,13500,135200,135525,135536),
            ),
            BinToPDFmap = cms.vstring("gaussPlusExpo")
        ))
        setattr(process.TnP_Trigger.Efficiencies, T+"_run_eta", cms.PSet(
            EfficiencyCategoryAndState = cms.vstring(T,"pass"),
            UnbinnedVariables = cms.vstring("mass"),
            BinnedVariables = cms.PSet(
                Glb = cms.vstring("true"),
                pt  = cms.vdouble( 3, 15),
                eta = cms.vdouble(-2.1,-1.0,1.0, 2.1),
                run = cms.vdouble(132440,133443,13500,135200,135525,135536),
            ),
            BinToPDFmap = cms.vstring("gaussPlusExpo")
        ))



process.p = cms.Path(
    process.TnP_Trigger
)


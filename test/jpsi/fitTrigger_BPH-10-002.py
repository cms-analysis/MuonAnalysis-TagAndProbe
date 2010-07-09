import FWCore.ParameterSet.Config as cms

### USAGE:
###    cmsRun fitTrigger_ICHEP.py <scenario>
### scenarios:
###   - data_all:    will fit tnpJPsi_Data.root with bins suitable for the current data
###   - datalike_mc: will fit tnpJPsi_{JPsiMuMu,ppMuX}_Spring10_0.117pb.root MC but
###                  with same config as data

import sys
args = sys.argv[1:]
if (sys.argv[0] == "cmsRun"): args =sys.argv[2:]
scenario = "data_all"
if len(args) > 0: scenario = args[0]
print "Will run scenario ", scenario 

process = cms.Process("TagProbe")

process.load('FWCore.MessageService.MessageLogger_cfi')

process.source = cms.Source("EmptySource")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

Template = cms.EDAnalyzer("TagProbeFitTreeAnalyzer",
    NumCPU = cms.uint32(1),
    SaveWorkspace = cms.bool(False),

    Variables = cms.PSet(
        mass   = cms.vstring("Tag-Probe Mass", "2.8", "3.5", "GeV/c^{2}"),
        pt     = cms.vstring("Probe p_{T}", "0", "1000", "GeV/c"),
        abseta = cms.vstring("Probe |#eta|", "0", "2.5", ""),
        #eta = cms.vstring("Probe #eta", "-2.5", "2.5", ""),
        tag_pt = cms.vstring("Tag p_{T}", "2.6", "1000", "GeV/c"),
    ),

    Categories = cms.PSet(
        # Efficiencies to measure
        L1DoubleMuOpen = cms.vstring("L1DoubleMuOpen", "dummy[pass=1,fail=0]"),
        Mu3 = cms.vstring("Mu3", "dummy[pass=1,fail=0]"),
        # Denominator
        Glb = cms.vstring("Glb", "dummy[pass=1,fail=0]"),
        TM = cms.vstring("TM", "dummy[pass=1,fail=0]"),
        POG_TMLSAT = cms.vstring("POG_TMLSAT", "dummy[pass=1,fail=0]"),
        # Other constraints
        tag_Mu3 = cms.vstring("tag_Mu3", "dummy[pass=1,fail=0]"),
        Acc_JPsi   = cms.vstring("tag_Mu3", "dummy[pass=1,fail=0]"),
        mcTrue = cms.vstring("MC true", "dummy[true=1,false=0]"),
    ),

    PDFs = cms.PSet(
        gaussPlusExpo = cms.vstring(
            "Gaussian::signal(mass, mean[3.1,3.0,3.2], sigma[0.05,0.02,0.1])",
            "Exponential::backgroundPass(mass, lp[0,-5,5])",
            "Exponential::backgroundFail(mass, lf[0,-5,5])",
            "efficiency[0.9,0,1]",
            "signalFractionInPassing[0.9]"
        )
    ),

    Efficiencies = cms.PSet(), # will be filled later
)

CONSTRAINTS = cms.PSet(
    tag_Mu3 = cms.vstring("pass"),
    Acc_JPsi = cms.vstring("pass"),
)
PT_ETA = cms.PSet(
    CONSTRAINTS,
    pt  = cms.vdouble( 0.8, 1.5, 2, 2.5, 3.3, 4, 5, 8, 30 ),
    abseta = cms.vdouble( 0,   1.3, 2.4 ),
    #eta = cms.vdouble( -2.4, -1.3, 0, 1.3, 2.4 ),
)

#PREFIX="/data/gpetrucc/7TeV/tnp/trees/dev-jul02/"
#PREFIX="/afs/cern.ch/user/g/gpetrucc/scratch0/tnp/trees-07.07.2010-10am/"
PREFIX="/data/gpetrucc/7TeV/tnp/trees/dev-jul08-v1/"
process.TnP_Trigger = Template.clone(
    InputFileNames = cms.vstring(
        PREFIX+'tnpJPsi_Data_run132440to135735.root',
        PREFIX+'tnpJPsi_Data_run136033to137028.root',
        PREFIX+'tnpJPsi_Data_run138560to138751.root',
        PREFIX+'tnpJPsi_Data_run138919to139100.root',
        PREFIX+'tnpJPsi_Data_run139102to139195.root',
        PREFIX+'tnpJPsi_Data_run139239to139365.root',
        PREFIX+'tnpJPsi_Data_run139368to139400.root',
        PREFIX+'tnpJPsi_Data_run139407to139459.root',
    ),
    InputDirectoryName = cms.string("histoTrigger"),
    InputTreeName = cms.string("fitter_tree"),
    OutputFileName = cms.string("TnP_BPH-10-002_Trigger_%s.root" % scenario),
    Efficiencies = cms.PSet(),
)

if scenario == "data_all":
    process.TnP_Trigger.binsForMassPlots = cms.uint32(35)

if scenario == "datalike_mc":
    process.TnP_Trigger.InputFileNames = [
        PREFIX+"tnpJPsi_MC_JPsiToMuMu_0.122pb.root",
        PREFIX+"tnpJPsi_MC_ppMuX_0.122pb.root",
    ]

if scenario == "signal_mc":
    process.TnP_Trigger.InputFileNames = [ PREFIX+"tnpJPsi_MC_JPsiToMuMu_1.0pb.root" ]


for T in [ "L1DoubleMuOpen", "Mu3" ]:
    #for M in ["Glb", "TM_Excl", "TM_Incl"]:
    for M in ["Glb", "TM_Incl"]:
        #for X,B in [("barrel",PT_BARREL), ("endcaps",PT_ENDCAPS)]:
            X="abseta";
            BINNEDVARS = PT_ETA.clone()
            if M == "Glb":     setattr(BINNEDVARS, "Glb", cms.vstring("pass"))
            if M == "TM_Incl": setattr(BINNEDVARS, "TM", cms.vstring("pass"))
            if M == "TM_Excl": 
                setattr(BINNEDVARS, "TM",  cms.vstring("pass"))
                setattr(BINNEDVARS, "Glb", cms.vstring("fail"))
            if scenario != "signal_all":
                setattr(process.TnP_Trigger.Efficiencies, M+"_To_"+T+"_pt_"+X, cms.PSet(
                    EfficiencyCategoryAndState = cms.vstring(T,"pass"),
                    UnbinnedVariables = cms.vstring("mass"),
                    BinnedVariables = BINNEDVARS,
                    BinToPDFmap = cms.vstring("gaussPlusExpo")
                ))
            if scenario == "datalike_mc" or scenario == "signal_all":
                setattr(process.TnP_Trigger.Efficiencies, M+"_To_"+T+"_pt_"+X+"_mcTrue", cms.PSet(
                    EfficiencyCategoryAndState = cms.vstring(T,"pass"),
                    UnbinnedVariables = cms.vstring("mass"),
                    BinnedVariables = BINNEDVARS.clone(mcTrue = cms.vstring("true"))
                ))

process.p = cms.Path(
    process.TnP_Trigger
)


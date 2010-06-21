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
        tag_pt = cms.vstring("Tag p_{T}", "2.6", "1000", "GeV/c"),
    ),

    Categories = cms.PSet(
        # Efficiencies to measure
        L1DoubleMuOpen = cms.vstring("L1DoubleMuOpen", "dummy[pass=1,fail=0]"),
        Mu3 = cms.vstring("Mu3", "dummy[pass=1,fail=0]"),
        # Denominator
        POG_Glb = cms.vstring("POG_Glb", "dummy[pass=1,fail=0]"),
        POG_GlbPT = cms.vstring("POG_GlbPT", "dummy[pass=1,fail=0]"),
        POG_TMA = cms.vstring("POG_TMA", "dummy[pass=1,fail=0]"),
        POG_TMLSAT = cms.vstring("POG_TMLSAT", "dummy[pass=1,fail=0]"),
        Cal = cms.vstring("Cal", "dummy[pass=1,fail=0]"),
        VBTFLike = cms.vstring("VBTFLike", "dummy[pass=1,fail=0]"),
        # Other constraints
        tag_Mu3 = cms.vstring("tag_Mu3", "dummy[pass=1,fail=0]"),
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
)
PT_ETA_BINS_DET = cms.PSet(
    CONSTRAINTS,
    pt     = cms.vdouble(  2.0, 3.0, 5.0, 15),
    abseta = cms.vdouble(  0.0, 0.8, 1.2, 1.6, 2.1, 2.4), 
)
PT_ETA_BINS_PH = cms.PSet(
    CONSTRAINTS,
    pt     = cms.vdouble(  2.0, 3.0, 5.0, 7.0, 15),
    abseta = cms.vdouble(  0.0, 1.2, 2.4),
)
PT_ETA_BINS = PT_ETA_BINS_DET




process.TnP_Trigger = Template.clone(
    InputFileNames = cms.vstring(
        #'/afs/cern.ch/user/g/gpetrucc/scratch0/tnp/tnpJPsi_Data_run132440to134987.root', # NO: taken online with 3.5.4
        '/afs/cern.ch/user/g/gpetrucc/scratch0/tnp/tnpJPsi_Data_run13509to135175.root',
        '/afs/cern.ch/user/g/gpetrucc/scratch0/tnp/tnpJPsi_Data_run135445to135575.root',
        '/afs/cern.ch/user/g/gpetrucc/scratch0/tnp/tnpJPsi_Data_run135735.root',
        '/afs/cern.ch/user/g/gpetrucc/scratch0/tnp/tnpJPsi_Data_run136033to136082.root',
        '/afs/cern.ch/user/g/gpetrucc/scratch0/tnp/tnpJPsi_Data_run136087to136119.root',
        '/afs/cern.ch/user/g/gpetrucc/scratch0/tnp/tnpJPsi_Data_run137027to137028.root',
    ),
    InputDirectoryName = cms.string("histoTrigger"),
    InputTreeName = cms.string("fitter_tree"),
    OutputFileName = cms.string("TnP_ICHEP_Trigger_%s.root" % scenario),
    Efficiencies = cms.PSet(),
)

if scenario == "data_all":
    process.TnP_Trigger.binsForMassPlots = cms.uint32(35)

if scenario == "datalike_mc":
    process.TnP_Trigger.InputFileNames = [
        "/afs/cern.ch/user/g/gpetrucc/scratch0/tnp/tnpJPsi_JPsiMuMu_Spring10_0.117pb.root",
        "/afs/cern.ch/user/g/gpetrucc/scratch0/tnp/tnpJPsi_ppMuX_Spring10_0.117pb.root"
    ]



for T in [ "L1DoubleMuOpen", "Mu3" ]:
    #for M in ["POG_Glb", "POG_GlbPT", "POG_TMA", "POG_TMLSAT"]:
    for M in ["POG_Glb", "VBTFLike", "Cal"]:
        for BN,BV in (('abseta',PT_ETA_BINS_DET),('pt',PT_ETA_BINS_PH)):
            BINNEDVARS = BV.clone()
            setattr(BINNEDVARS, M, cms.vstring("pass"))
            setattr(process.TnP_Trigger.Efficiencies, M+"_To_"+T+"_"+BN, cms.PSet(
                EfficiencyCategoryAndState = cms.vstring(T,"pass"),
                UnbinnedVariables = cms.vstring("mass"),
                BinnedVariables = BINNEDVARS,
                BinToPDFmap = cms.vstring("gaussPlusExpo")
            ))
            if scenario == "datalike_mc":
                setattr(process.TnP_Trigger.Efficiencies, M+"_To_"+T+"_"+BN+"_mcTrue", cms.PSet(
                    EfficiencyCategoryAndState = cms.vstring(T,"pass"),
                    UnbinnedVariables = cms.vstring("mass"),
                    BinnedVariables = BINNEDVARS.clone(mcTrue = cms.vstring("true"))
                ))
            if T != "L1DoubleMuOpen":
                setattr(process.TnP_Trigger.Efficiencies, M+"_To_"+T+"overL1_"+BN, cms.PSet(
                    EfficiencyCategoryAndState = cms.vstring(T,"pass"),
                    UnbinnedVariables = cms.vstring("mass"),
                    BinnedVariables = BINNEDVARS.clone(L1DoubleMuOpen = cms.vstring("pass")),
                    BinToPDFmap = cms.vstring("gaussPlusExpo")
                ))
                if scenario == "datalike_mc":
                    setattr(process.TnP_Trigger.Efficiencies, M+"_To_"+T+"overL1_"+BN+"_mcTrue", cms.PSet(
                        EfficiencyCategoryAndState = cms.vstring(T,"pass"),
                        UnbinnedVariables = cms.vstring("mass"),
                        BinnedVariables = PT_ETA_BINS.clone(L1DoubleMuOpen = cms.vstring("pass"), mcTrue = cms.vstring("true"))
                    )) 

process.p = cms.Path(
    process.TnP_Trigger
)


import FWCore.ParameterSet.Config as cms

### USAGE:
###    cmsRun fitTrigger_Z.py <scenario>
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

process.TnP_Trigger = cms.EDAnalyzer("TagProbeFitTreeAnalyzer",
    NumCPU = cms.uint32(1),
    SaveWorkspace = cms.bool(False),

    InputFileNames = cms.vstring('tnpZ_Data.root'),
    InputTreeName = cms.string("fitter_tree"),
    InputDirectoryName = cms.string("tpTree"),
    OutputFileName = cms.string("TnP_Z_Trigger_%s.root" % scenario),

    Variables = cms.PSet(
        mass = cms.vstring("Tag-Probe Mass", "50", "130", "GeV/c^{2}"),
        pt     = cms.vstring("Probe p_{T}", "0", "1000", "GeV/c"),
        eta    = cms.vstring("Probe |#eta|", "-2.5", "2.5", ""),
        abseta = cms.vstring("Probe |#eta|", "0", "2.5", ""),
        tag_pt = cms.vstring("Tag p_{T}", "2.6", "1000", "GeV/c"),
    ),

    Categories = cms.PSet(
        Calo = cms.vstring("POG_Glb",  "dummy[pass=1,fail=0]"),
        Glb  = cms.vstring("POG_Glb",  "dummy[pass=1,fail=0]"),
        VBTF = cms.vstring("VBTFLike", "dummy[pass=1,fail=0]"),
        Isol = cms.vstring("MC true",  "dummy[pass=1,fail=0]"),
        Mu9  = cms.vstring("MC true",  "dummy[pass=1,fail=0]"),
    ),

    PDFs = cms.PSet(
        gaussPlusExpo = cms.vstring(
            "Voigtian::signal(mass, mean[90,80,100], width[2.495], sigma[3,1,20])",
            "Exponential::backgroundPass(mass, lp[0,-5,5])",
            "Exponential::backgroundFail(mass, lf[0,-5,5])",
            "efficiency[0.9,0,1]",
            "signalFractionInPassing[0.9]"
        )
    ),

    Efficiencies = cms.PSet(), # will be filled later
)
ONE_BIN = cms.PSet(
    pt     = cms.vdouble(15, 100),
    abseta = cms.vdouble( 0, 2.1)
)
PT_BINS  = ONE_BIN.clone(pt = cms.vdouble(15, 35, 100))
ETA_BINS = ONE_BIN.clone(abseta = cms.vdouble(0, 1.2, 2.1))


if scenario == "data_all":
    process.TnP_Trigger.binsForMassPlots = cms.uint32(20)

if scenario == "datalike_mc":
    process.TnP_Trigger.InputFileNames = [ "tnpZ_MC.root", ]


ALLBINS=[("all",ONE_BIN),("pt",PT_BINS),("eta",ETA_BINS)]
for (T,M) in [ ("Mu9","Track"), ("Mu9","Glb"), ("Mu9", "VBTF"), ("Mu9", "VBTF_Isol")]:
        for BN,BV in ALLBINS:
            BINNEDVARS = BV.clone()
            if M == "VBTF_Isol":
                setattr(BINNEDVARS, "VBTF", cms.vstring("pass"))
                setattr(BINNEDVARS, "Isol", cms.vstring("pass"))
            elif M != "Track": 
                setattr(BINNEDVARS, M, cms.vstring("pass"))
            setattr(process.TnP_Trigger.Efficiencies, M+"_To_"+T+"_"+BN, cms.PSet(
                EfficiencyCategoryAndState = cms.vstring(T,"pass"),
                UnbinnedVariables = cms.vstring("mass"),
                BinnedVariables = BINNEDVARS,
                BinToPDFmap = cms.vstring("gaussPlusExpo")
            ))
for  X,B in ALLBINS:
    setattr(process.TnP_Trigger.Efficiencies, "Track_To_VBTF_Mu9_"+X, cms.PSet(
        EfficiencyCategoryAndState = cms.vstring("VBTF","pass","Mu9","pass"),
        UnbinnedVariables = cms.vstring("mass"),
        BinnedVariables = B,
        BinToPDFmap = cms.vstring("gaussPlusExpo")
    ))
    setattr(process.TnP_Trigger.Efficiencies, "Track_To_VBTF_Isol_Mu9_"+X, cms.PSet(
        EfficiencyCategoryAndState = cms.vstring("VBTF","pass","Isol","pass","Mu9","pass"),
        UnbinnedVariables = cms.vstring("mass"),
        BinnedVariables = B,
        BinToPDFmap = cms.vstring("gaussPlusExpo")
    ))


process.p = cms.Path(
    process.TnP_Trigger
)


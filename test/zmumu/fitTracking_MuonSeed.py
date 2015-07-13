import FWCore.ParameterSet.Config as cms

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

CONSTRAINTS = cms.PSet(
    #staValidStations = cms.vdouble(0.5,99),
)

ONE_BIN = cms.PSet(CONSTRAINTS,
    pt = cms.vdouble( 15, 120 ),
    eta = cms.vdouble(-2.4, 2.4),
)
ONE_BIN1 = cms.PSet(CONSTRAINTS,
    pt = cms.vdouble( 15, 120 ),
    eta = cms.vdouble(-2.1, 2.1),
)
TWO_BINS = cms.PSet(CONSTRAINTS,
    pt = cms.vdouble( 15, 120 ),
    abseta = cms.vdouble(0, 1.2, 2.4),
)

ETA_BINS = cms.PSet(CONSTRAINTS,
    pt = cms.vdouble( 15, 120 ),
    eta = cms.vdouble(*[-2.4+0.2*x for x in range(0,25)])
)
AETA_BINS = cms.PSet(CONSTRAINTS,
    pt = cms.vdouble( 15, 120 ),
    abseta = cms.vdouble(0, 0.6, 1.1, 1.6, 2.1, 2.4),
)
PT_AETA_BINS1 = cms.PSet(CONSTRAINTS,
    pt = cms.vdouble( 5, 8, 10, 12, 15, 20, 120 ),
    abseta = cms.vdouble(0, 1.2, 2.5),
)


process.TnP_Muon_Seed = cms.EDAnalyzer("TagProbeFitTreeAnalyzer",
    ## Input, output 
    InputFileNames = cms.vstring(), ## can put more than one
    OutputFileName = cms.string("TnP_Muon_Seed_%s.root" % scenario),
    InputTreeName = cms.string("fitter_tree"), 
    InputDirectoryName = cms.string("tpTree"),  
    ## Variables for binning
    Variables = cms.PSet(
        mass   = cms.vstring("Tag-muon Mass", "76", "125", "GeV/c^{2}"),
        pt     = cms.vstring("muon p_{T}", "0", "1000", "GeV/c"),
        eta    = cms.vstring("Probe #eta", "-2.5", "2.5", ""),
        abseta = cms.vstring("muon |#eta|", "0", "2.5", ""),
        staValidStations = cms.vstring("Valid stations in muon system", "-2", "10", "cm"),
    ),
    ## Flags you want to use to define numerator and possibly denominator
    Categories = cms.PSet(
        MuIDForOutsideInTk = cms.vstring("OITK Seed", "dummy[pass=1,fail=0]"),
        #tag_IsoMu24_eta2p1 = cms.vstring("Trigger", "dummy[pass=1,fail=0]"),
    ),
    ## What to fit
    Efficiencies = cms.PSet(),
    ## PDF for signal and background (double voigtian + exponential background)
    PDFs = cms.PSet(
        voigtPlusExpo = cms.vstring(
            "Voigtian::signal(mass, meanP[90,80,100], width[2.495], sigma[5,1,12])",
            "Exponential::backgroundPass(mass, lp[0,-5,5])",
            "Exponential::backgroundFail(mass, lf[0,-5,5])",
            "efficiency[0.9,0,1]",
            "signalFractionInPassing[0.9]"
        ),
        vpvPlusExpo = cms.vstring(
            "Voigtian::signal1(mass, mean1[90,80,100], width[2.495], sigma1[2,1,3])",
            "Voigtian::signal2(mass, mean2[90,80,100], width,        sigma2[4,2,10])",
            "SUM::signal(vFrac[0.8,0,1]*signal1, signal2)",
            "Exponential::backgroundPass(mass, lp[-0.1,-1,0.1])",
            "Exponential::backgroundFail(mass, lf[-0.1,-1,0.1])",
            "efficiency[0.9,0,1]",
            "signalFractionInPassing[0.9]"
        ),
    ),
    ## How to do the fit
    binnedFit = cms.bool(True),
    binsForFit = cms.uint32(40),
    saveDistributionsPlot = cms.bool(False),
    NumCPU = cms.uint32(1), ## leave to 1 for now, RooFit gives funny results otherwise
    SaveWorkspace = cms.bool(False),
)

if "TestMC_" in scenario:
    process.TnP_Muon_Seed.InputFileNames = [ "tnpZ_MC_%s.root" % ("_".join(scenario.replace("TestMC_","").split("_")[:-1]),) ]

common = cms.PSet(
    EfficiencyCategoryAndState = cms.vstring("MuIDForOutsideInTk","pass"),
    UnbinnedVariables = cms.vstring("mass"),
    BinToPDFmap = cms.vstring("voigtPlusExpo")
)

BINS =  [ ]
BINS += [ ("avg", "", ONE_BIN ) ]
BINS += [ ("av1", "", ONE_BIN1 ) ]
BINS += [ ("two", "_two", TWO_BINS ) ]
BINS += [ ("_eta", "_eta", ETA_BINS ) ]
BINS += [ ("aeta", "_aeta", AETA_BINS ) ]
BINS += [ ("pa1", "_pa1", PT_AETA_BINS1 ) ]
for (n,post,B) in BINS:
    if n in scenario:
        B0 = B.clone()
        B1 = B.clone(staValidStations = cms.vdouble(0.5,99))
        setattr(process.TnP_Muon_Seed.Efficiencies, "effS"+post,  cms.PSet(common, BinnedVariables = B0))
        setattr(process.TnP_Muon_Seed.Efficiencies, "effSP"+post, cms.PSet(common, BinnedVariables = B1))

process.p = cms.Path(process.TnP_Muon_Seed)

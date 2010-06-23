import FWCore.ParameterSet.Config as cms

import sys
args = sys.argv[1:]
if (sys.argv[0] == "cmsRun"): args =sys.argv[2:]
scenario = "data_all"
if len(args) > 0: scenario = args[0]
print "Will run scenario ", scenario 


CONSTRAINTS = cms.PSet(
    hasValidHits = cms.vstring("pass"),
    #tag_Mu3 = cms.vstring("pass"),
    tag_L1DoubleMuOpen = cms.vstring("pass"),
)
ONE_BIN = cms.PSet(CONSTRAINTS,
    pt = cms.vdouble( 0, 15 ),
    eta = cms.vdouble(-2.5, 2.5),
)
ETA_BINS = ONE_BIN.clone(
    eta = cms.vdouble(-2.5, -1.1, 1.1, 2.5)
)
ETA_PHI_BINS =  ETA_BINS.clone(
   phi = cms.vdouble(*[3.1416*i/3.0 for i in range(-3,4)]), # 6 bins
)

process = cms.Process("TagProbe")

process.load('FWCore.MessageService.MessageLogger_cfi')

process.source = cms.Source("EmptySource")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

Template = cms.EDAnalyzer("TagProbeFitTreeAnalyzer",
    NumCPU = cms.uint32(1),
    SaveWorkspace = cms.bool(True),

    Variables = cms.PSet(
        mass   = cms.vstring("Tag-Probe Mass", "2.0", "4.3", "GeV/c^{2}"),
        pt     = cms.vstring("Probe p_{T}", "0", "1000", "GeV/c"),
        p      = cms.vstring("Probe p", "0", "1000", "GeV/c"),
        eta    = cms.vstring("Probe #eta", "-2.5", "2.5", ""),
        abseta = cms.vstring("Probe |#eta|", "0", "2.5", ""),
        phi    = cms.vstring("Probe #phi", "-3.1416", "3.1416", ""),
        tag_pt = cms.vstring("Tag p_{T}", "2.6", "1000", "GeV/c"),
        match_deltaR   = cms.vstring("Match #Delta R",    "0", "1000", ""),
        match_deltaEta = cms.vstring("Match #Delta #eta", "0", "1000", ""),
        match_deltaR_NoJPsi   = cms.vstring("Unmatch #Delta R",    "0", "1000", ""),
        match_deltaEta_NoJPsi = cms.vstring("Unmatch #Delta #eta", "0", "1000", ""),
        match_deltaR_NoBestJPsi   = cms.vstring("Unmatch #Delta R",    "0", "1000", ""),
        match_deltaEta_NoBestJPsi = cms.vstring("Unmatch #Delta #eta", "0", "1000", ""),
    ),

    Categories = cms.PSet(
        passing            = cms.vstring("passing",           "dummy[pass=1,fail=0]"),
        passingNoJPsi      = cms.vstring("passingNoJPsi",     "dummy[pass=1,fail=0]"),
        passingNoBestJPsi  = cms.vstring("passingNoBestJPsi", "dummy[pass=1,fail=0]"),
        ## Quality cut
        hasValidHits = cms.vstring("hasValidHits",  "dummy[pass=1,fail=0]"),
        ## Trigger
        tag_Mu3 = cms.vstring("tag_HLTMu3", "dummy[pass=1,fail=0]"),
        tag_L1DoubleMuOpen = cms.vstring("tag_HLTL1DoubleMuOpen", "dummy[pass=1,fail=0]"),
        ## MC
        mcTrue = cms.vstring("MC true", "dummy[true=1,false=0]"),
    ),

    Expressions = cms.PSet(),
    Cuts = cms.PSet(),

    PDFs = cms.PSet(
        gaussPlusCubic = cms.vstring(
            "Gaussian::signal(mass, mean[3.1,3.0,3.2], sigma[0.2,0.05,0.5])",
            "Chebychev::backgroundPass(mass, {c1[0,-1,1], c2[0,-1,1], c3[0,-1,1]})",
            "Chebychev::backgroundFail(mass, {c1,c2,c3})",
            "efficiency[0.9,0,1]",
            "signalFractionInPassing[0.5]"
        ),
        gaussPlusFloatCubic = cms.vstring(
            "Gaussian::signal(mass, mean[3.1,3.0,3.2], sigma[0.2,0.05,0.5])",
            "Chebychev::backgroundPass(mass, {c1p[0,-1,1], c2p[0,-1,1], c3p[0,-1,1]})",
            "Chebychev::backgroundFail(mass, {c1f[0,-1,1], c2f[0,-1,1], c3f[0,-1,1]})",
            "efficiency[0.9,0,1]",
            "signalFractionInPassing[0.5]"
        )

    )
)

tofit = []
#for (dr,de) in [ (1.0,0.4), (0.7,0.3), (0.5,0.2), (0.4,0.2), (0.3,0.15), (0.2,0.1), (0.1,0.1), (0.05,0.05) ): 
for (dr,de) in [ (0.5,0.2) ]: 
    tofit.append("dr%03de%03d" % (100*dr, 100*de))
    #for X in ("", "_NoJPsi","_NoBestJPsi"):
    for X in ("", "_NoBestJPsi"):
        name = "dr%03de%03d%s" % (100*dr, 100*de, X.replace("_",""))
        cut  = "match_deltaR%s < %f && match_deltaEta%s < %f" % (X, dr, X, de)
        setattr(Template.Expressions, name, cms.vstring(cut,cut,"match_deltaR"+X, "match_deltaEta"+X))
for N in Template.Expressions.parameterNames_(): 
    setattr(Template.Cuts, "cut_"+N, cms.vstring(N, N, "0.5"))

process.TnP_Tracking = Template.clone(
    InputFileNames = cms.vstring("tnpJPsi_JPsiMuMu_Spring10_0.1pb.root"),
    InputDirectoryName = cms.string("histoTracking"),
    InputTreeName = cms.string("fitter_tree"),
    OutputFileName = cms.string("TnP_Tracking_%s.root" % scenario),
    Efficiencies = cms.PSet()
)

if scenario == "data_all":
    process.TnP_Tracking.InputFileNames = cms.vstring(
        '/afs/cern.ch/user/g/gpetrucc/scratch0/tnp/tnpJPsi_Data_run132440to134987.root',
        '/afs/cern.ch/user/g/gpetrucc/scratch0/tnp/tnpJPsi_Data_run13509to135175.root',
        '/afs/cern.ch/user/g/gpetrucc/scratch0/tnp/tnpJPsi_Data_run135445to135575.root',
        '/afs/cern.ch/user/g/gpetrucc/scratch0/tnp/tnpJPsi_Data_run135735.root',
        '/afs/cern.ch/user/g/gpetrucc/scratch0/tnp/tnpJPsi_Data_run136033to136082.root',
        '/afs/cern.ch/user/g/gpetrucc/scratch0/tnp/tnpJPsi_Data_run136087to136119.root',
        '/afs/cern.ch/user/g/gpetrucc/scratch0/tnp/tnpJPsi_Data_run137027to137028.root',
    )
    process.TnP_Tracking.binsForMassPlots = cms.uint32(23)
elif scenario == "datalike_mc":
    process.TnP_Tracking.InputFileNames = [
        "/afs/cern.ch/user/g/gpetrucc/scratch0/tnp/tnpJPsi_JPsiMuMu_Spring10_0.117pb.root",
        "/afs/cern.ch/user/g/gpetrucc/scratch0/tnp/tnpJPsi_ppMuX_Spring10_0.117pb.root"
    ]

sampleToPdfMap = { "": "gaussPlusCubic", "NoJPsi":"gaussPlusFloatCubic", "NoBestJPsi":"gaussPlusFloatCubic" }
for M in tofit:
    #for X in ("","NoJPsi","NoBestJPsi"):
    for X in ("","NoBestJPsi"):
        common = cms.PSet(
            EfficiencyCategoryAndState = cms.vstring("cut_"+M+X,"above"),
            UnbinnedVariables = cms.vstring("mass"),
            BinToPDFmap = cms.vstring(sampleToPdfMap[X])
        )
        setattr(process.TnP_Tracking.Efficiencies, "eff_"    +M+X, cms.PSet(common, BinnedVariables = ONE_BIN))
        if True:
            setattr(process.TnP_Tracking.Efficiencies, "eff_eta_"+M+X, cms.PSet(common, BinnedVariables = ETA_BINS))

process.p = cms.Path(
    process.TnP_Tracking
)


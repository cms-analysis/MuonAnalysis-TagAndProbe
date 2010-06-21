import FWCore.ParameterSet.Config as cms

### USAGE:
###    cmsRun fitMuonID_ICHEP.py <scenario>
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
    SaveWorkspace = cms.bool(True),

    Variables = cms.PSet(
        mass = cms.vstring("Tag-Probe Mass", "2.8", "3.5", "GeV/c^{2}"),
        pt = cms.vstring("Probe p_{T}", "0", "1000", "GeV/c"),
        #p = cms.vstring("Probe p", "0", "1000", "GeV/c"),
        abseta = cms.vstring("Probe |#eta|", "0", "2.5", ""),
        #phi = cms.vstring("Probe #phi", "-3.1416", "3.1416", ""),
        tag_pt = cms.vstring("Tag p_{T}", "2.6", "1000", "GeV/c"),
    ),

    Categories = cms.PSet(
        POG_TMLSAT = cms.vstring("POG_TMLSAT", "dummy[pass=1,fail=0]"),
        VBTFLike = cms.vstring("VBTFLike", "dummy[pass=1,fail=0]"),
        POG_Glb = cms.vstring("POG_Glb", "dummy[pass=1,fail=0]"),
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
PT_ETA_BINS = cms.PSet(
    CONSTRAINTS,
    pt     = cms.vdouble(  0.5, 2.0, 3.0, 5.0, 7.0, 15.0),
    abseta = cms.vdouble(  0.0, 1.2, 2.4)
)


process.TnP_MuonID = Template.clone(
    InputFileNames = cms.vstring(
        '/afs/cern.ch/user/g/gpetrucc/scratch0/tnp/tnpJPsi_Data_run132440to134987.root',
        '/afs/cern.ch/user/g/gpetrucc/scratch0/tnp/tnpJPsi_Data_run13509to135175.root',
        '/afs/cern.ch/user/g/gpetrucc/scratch0/tnp/tnpJPsi_Data_run135445to135575.root',
        '/afs/cern.ch/user/g/gpetrucc/scratch0/tnp/tnpJPsi_Data_run135735.root',
        '/afs/cern.ch/user/g/gpetrucc/scratch0/tnp/tnpJPsi_Data_run136033to136082.root',
        '/afs/cern.ch/user/g/gpetrucc/scratch0/tnp/tnpJPsi_Data_run136087to136119.root',
        '/afs/cern.ch/user/g/gpetrucc/scratch0/tnp/tnpJPsi_Data_run137027to137028.root',
    ),
    InputTreeName = cms.string("fitter_tree"),
    InputDirectoryName = cms.string("histoMuFromCal"),
    OutputFileName = cms.string("TnP_ICHEP_MuonID_%s.root" % scenario),
    Efficiencies = cms.PSet(),
)

if scenario == "data_all":
    process.TnP_MuonID.binsForMassPlots = cms.uint32(35)

if scenario == "datalike_mc":
    process.TnP_MuonID.InputFileNames = [
        "/afs/cern.ch/user/g/gpetrucc/scratch0/tnp/tnpJPsi_JPsiMuMu_Spring10_0.117pb.root",
        "/afs/cern.ch/user/g/gpetrucc/scratch0/tnp/tnpJPsi_ppMuX_Spring10_0.117pb.root"
    ]



for T in [ "POG_Glb", "POG_TMLSAT", "VBTFLike" ]:
    setattr(process.TnP_MuonID.Efficiencies, T+"_pt_abseta", cms.PSet(
        EfficiencyCategoryAndState = cms.vstring(T,"pass"),
        UnbinnedVariables = cms.vstring("mass"),
        BinnedVariables = PT_ETA_BINS,
        BinToPDFmap = cms.vstring("gaussPlusExpo")
    ))
    if scenario == "datalike_mc":
        setattr(process.TnP_MuonID.Efficiencies, T+"_pt_abseta_mcTrue", cms.PSet(
            EfficiencyCategoryAndState = cms.vstring(T,"pass"),
            UnbinnedVariables = cms.vstring("mass"),
            BinnedVariables = PT_ETA_BINS.clone(mcTrue = cms.vstring("true"))
        )) 

process.TnP_MuonID_Tk = process.TnP_MuonID.clone(
    InputDirectoryName = cms.string("histoMuFromTk"),
    OutputFileName = cms.string("TnP_ICHEP_MuonID_FromTK_%s.root" % scenario),
)

process.p = cms.Path(
    process.TnP_MuonID    +
    process.TnP_MuonID_Tk
)


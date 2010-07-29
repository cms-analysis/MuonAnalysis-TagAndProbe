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
doTk, doCal = True, True
if scenario.find("_tk") != -1: 
    print "Running only eff from tracks"
    doCal = False
    scenario = scenario.replace("_tk","")
if scenario.find("_cal") != -1: 
    print "Running only eff from calo"
    doTk = False
    scenario = scenario.replace("_cal","")


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
        #eta = cms.vstring("Probe #eta", "-2.5", "2.5", ""),
        abseta = cms.vstring("Probe |#eta|", "0", "2.5", ""),
        #phi = cms.vstring("Probe #phi", "-3.1416", "3.1416", ""),
        tag_pt = cms.vstring("Tag p_{T}", "2.6", "1000", "GeV/c"),
        pair_Nvertices = cms.vstring("Number of vertices", "0", "999", ""),
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
    pt     = cms.vdouble(  0.5, 2.0, 2.5, 3.0, 4.0, 5.0, 7.0, 14.0),
    abseta = cms.vdouble(  0.0, 1.2, 2.4)
)
VTX_BINS = cms.PSet(
    CONSTRAINTS,
    pt     = cms.vdouble(  0.5, 2.0, 4.5, 14.0),
    abseta = cms.vdouble(  0.0, 1.2, 2.4),
    pair_Nvertices = cms.vdouble(0.5,1.5,2.5,3.5,4.5)
)



PREFIX="/afs/cern.ch/user/g/gpetrucc/scratch0/tnp/trees-14.07.2010/"
process.TnP_MuonID = Template.clone(
    InputFileNames = cms.vstring(
        PREFIX+'tnpJPsi_Data_run132440to135735.root',
        PREFIX+'tnpJPsi_Data_run136033to137028.root',
        PREFIX+'tnpJPsi_Data_run138560to138751.root',
        PREFIX+'tnpJPsi_Data_run138919to139100.root',
        PREFIX+'tnpJPsi_Data_run139102to139195.root',
        PREFIX+'tnpJPsi_Data_run139239to139365.root',
        PREFIX+'tnpJPsi_Data_run139368to139400.root',
        PREFIX+'tnpJPsi_Data_run139407to139459.root',
        PREFIX+'tnpJPsi_Data_run139779to139790.root',
    ),
    InputTreeName = cms.string("fitter_tree"),
    InputDirectoryName = cms.string("histoMuFromCal"),
    OutputFileName = cms.string("TnP_ICHEP_MuonID_%s.root" % scenario),
    Efficiencies = cms.PSet(),
)

if scenario == "data_all":
    process.TnP_MuonID.binsForMassPlots = cms.uint32(70)

if scenario == "datalike_mc":
    process.TnP_MuonID.InputFileNames = [
        PREFIX+"tnpJPsi_MC_JPsiToMuMu_0.122pb.root",
        PREFIX+"tnpJPsi_MC_ppMuX_0.122pb.root",
    ]

if scenario == "signal_mc":
    process.TnP_MuonID.InputFileNames = [ PREFIX+"tnpJPsi_MC_JPsiToMuMu_1.0pb.root" ]


ALLBINS=[("pt_abseta",PT_ETA_BINS)]
if scenario == "data_all": ALLBINS += [ ("vtx",VTX_BINS)]
for T in [ "POG_Glb", "POG_TMLSAT", "VBTFLike" ]:
    for  X,B in ALLBINS:
        setattr(process.TnP_MuonID.Efficiencies, T+"_"+X, cms.PSet(
            EfficiencyCategoryAndState = cms.vstring(T,"pass"),
            UnbinnedVariables = cms.vstring("mass"),
            BinnedVariables = B,
            BinToPDFmap = cms.vstring("gaussPlusExpo")
        ))
        if scenario == "datalike_mc" or scenario == "signal_mc":
            setattr(process.TnP_MuonID.Efficiencies, T+"_"+X+"_mcTrue", cms.PSet(
                EfficiencyCategoryAndState = cms.vstring(T,"pass"),
                UnbinnedVariables = cms.vstring("mass"),
                BinnedVariables = B.clone(mcTrue = cms.vstring("true"))
            )) 

process.TnP_MuonID_Tk = process.TnP_MuonID.clone(
    InputDirectoryName = cms.string("histoMuFromTk"),
    OutputFileName = cms.string("TnP_ICHEP_MuonID_FromTK_%s.root" % scenario),
)

#if doTk:  process.pTk  = cms.Path(process.TnP_MuonID_Tk)
if doCal: process.pCal = cms.Path(process.TnP_MuonID)


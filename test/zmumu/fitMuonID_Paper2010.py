import FWCore.ParameterSet.Config as cms
### USAGE:
###    cmsRun fitMuonID_Paper2010.py <scenario> [ <id> [ <binning1> ... <binningN> ] ]
###
### scenarios:
###   - data_all (default)  
###   - signal_mc

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
    saveDistributionsPlot = cms.bool(False),

    Variables = cms.PSet(
        mass = cms.vstring("Tag-muon Mass", "70", "130", "GeV/c^{2}"),
        pt = cms.vstring("muon p_{T}", "0", "1000", "GeV/c"),
        eta    = cms.vstring("muon #eta", "-2.5", "2.5", ""),
        abseta = cms.vstring("muon |#eta|", "0", "2.5", ""),
        phi    = cms.vstring("muon #phi at vertex", "-3.1416", "3.1416", ""),
        charge = cms.vstring("muon charge", "-2.5", "2.5", ""),
        tag_pt = cms.vstring("Tag p_{T}", "0", "1000", "GeV/c"),
        tag_nVertices = cms.vstring("Number of vertices", "0", "999", ""),
        isoTrk03Abs    = cms.vstring("Probe abs trk iso", "-2", "9999999", ""),
        tag_combRelIso = cms.vstring("Tag comb rel iso", "-2", "9999999", ""),
    ),

    Categories = cms.PSet(
        Glb   = cms.vstring("Global", "dummy[pass=1,fail=0]"),
        VBTF     = cms.vstring("VBTFLike", "dummy[pass=1,fail=0]"),
        VBTFold  = cms.vstring("VBTFLike", "dummy[pass=1,fail=0]"),
        TMOSL = cms.vstring("TMOneStationLoose", "dummy[pass=1,fail=0]"),
        TMOST = cms.vstring("TMOneStationTight", "dummy[pass=1,fail=0]"),
        PF    = cms.vstring("PF Muon", "dummy[pass=1,fail=0]"),
        TM    = cms.vstring("Tracker Muon", "dummy[pass=1,fail=0]"),
        TMA   = cms.vstring("Arbitrated Tracker Muon", "dummy[pass=1,fail=0]"),
        IsolTk3 = cms.vstring("Tk Abs Iso < 3", "dummy[pass=1,fail=0]"),
        Mu15  = cms.vstring("MC true", "dummy[true=1,false=0]"),
        mcTrue = cms.vstring("MC true", "dummy[true=1,false=0]"),
    ),

    PDFs = cms.PSet(
        voigtPlusExpo = cms.vstring(
            "Voigtian::signal(mass, mean[90,80,100], width[2.495], sigma[3,1,20])",
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
        )
    ),

    binnedFit = cms.bool(True),
    binsForFit = cms.uint32(40),

    Efficiencies = cms.PSet(), # will be filled later
)

PT_ETA_BINS = cms.PSet(
    pt     = cms.vdouble(  10, 20, 30, 40, 60, 100 ),
    abseta = cms.vdouble(  0.0, 1.2, 2.4)
)
ETA_BINS = cms.PSet(
    pt  = cms.vdouble(20,100),
    eta = cms.vdouble(-2.4, -2.1, -1.6, -1.1, -0.6, 0, 0.6, 1.1, 1.6, 2.1, 2.4),
)
ETA_PHI_BINS = ETA_BINS.clone(
   eta = cms.vdouble(-2.4, -1.6, -8.0, 0, 0.8, 1.6, 2.4),
   phi = cms.vdouble(*[3.1416*i/3.0 for i in range(-3,4)]), 
)
VTX_BINS  = cms.PSet(
    pt     = cms.vdouble(  20, 120 ),
    abseta = cms.vdouble(  0.0, 2.4),
    tag_nVertices = cms.vdouble(0.5,1.5,2.5,3.5,4.5,5.5,6.5)
)

ETA_BINS_FINE = cms.PSet(
    pt  = cms.vdouble(20,100),
    eta = cms.vdouble(*[x/10. for x in xrange(-24,25,1)]),
)
OVERALL = cms.PSet(
    pt  = cms.vdouble(20,100),
    abseta = cms.vdouble(0.0, 2.4),
)
OVERALL_ABSETA = cms.PSet(
    pt  = cms.vdouble(20,100),
    abseta = cms.vdouble(0.0, 1.2, 2.4),
)

OVERALL_ENDCAPS21 = cms.PSet(
    pt  = cms.vdouble(20,100),
    abseta = cms.vdouble(1.2, 2.1),
)


CHARGE = cms.PSet(
    pt     = cms.vdouble(20,100),
    abseta = cms.vdouble(0.0, 2.4),
    charge = cms.vdouble(-2,0,2),
)



#PREFIX="/data/gpetrucc/7TeV/tnp/2011.02.17/"
PREFIX="/data/gpetrucc/7TeV/tnp/2011-04/"
process.TnP_MuonID = Template.clone(
    InputFileNames = cms.vstring(PREFIX+'tnpZ_Data_Nov4B.root'),
    InputTreeName = cms.string("fitter_tree"),
    InputDirectoryName = cms.string("tpTree"),
    OutputFileName = cms.string("TnP_Paper2010_MuonID_%s.root" % scenario),
    Efficiencies = cms.PSet(),
)

if "39X" in scenario and "data" in scenario:
    process.TnP_MuonID.InputFileNames = [ PREFIX+"tnpZ_Data_Dec22B.root" ]
if scenario.find("signal_mc") != -1:
    process.TnP_MuonID.InputFileNames = [ PREFIX+"tnpZ_MC_DYPoweg.root" ]
if scenario.find("some_mc") != -1:
    process.TnP_MuonID.InputFileNames = [ PREFIX+"tnpZ_MCDYPoweg_38X_dzIso.334pb.root" ]
if scenario.find("39X_mc") != -1:
    process.TnP_MuonID.InputFileNames = [ PREFIX+"tnpZ_MC_39XDY_1.root" ]
if scenario.find("some_mc2") != -1:
    process.TnP_MuonID.InputFileNames = [ PREFIX+"tnpZ_DYToMMpowhegZ2_Fall10PU.root" ]
if scenario.find("39X_mc2") != -1:
    process.TnP_MuonID.InputFileNames = [ PREFIX+"tnpZ_DYToMMpowhegZ2_Winter10PU.root" ]
if scenario.find("realistic_mc") != -1:
    process.TnP_MuonID.InputFileNames = [ PREFIX+"tnpZ_MCDYPoweg_38X_dzIso.334pb.root", 
                                          PREFIX+"tnpZ_MC_WJetsPU.334pb.root", 
                                          PREFIX+"tnpZ_MC_QCDMuPt15PU.334pb.root", ]
if "tag35" in scenario:
    process.TnP_MuonID.Variables.tag_pt[1]='35'

IDS = [ "TMOST", "VBTF", "PF" ]
IDS += [ "Glb" ]
ALLBINS = [("pt_abseta",PT_ETA_BINS),("eta", ETA_BINS)]
ALLBINS += [ ("vtx",VTX_BINS)]
#ALLBINS+=[("eta_fine",ETA_BINS_FINE)]
ALLBINS += [("overall",OVERALL), ("charge",CHARGE), ("overall_abseta",OVERALL_ABSETA),("overall_endcaps21",OVERALL_ENDCAPS21)]
ALLBINS+=[("eta_phi",ETA_PHI_BINS)]

if len(args) > 1 and args[1] not in IDS: IDS += [ args[1] ]
for ID in IDS:
    if len(args) > 1 and ID != args[1]: continue
    for X,B in ALLBINS:
        if len(args) > 2 and X not in args[2:]: continue
        module = process.TnP_MuonID.clone(OutputFileName = cms.string("TnP_Paper2010_MuonID_%s_%s_%s.root" % (scenario, ID, X)))
        shape = "vpvPlusExpo"
        if X.find("eta") != -1 and X.find("abseta") == -1: shape = "voigtPlusExpo"
        if X.find("pt_abseta") != -1: module.Variables.mass[1]="77";
        if X.find("overall") != -1: module.binsForFit = 120
        DEN = B.clone(); num = ID;
        if "_from_" in ID:
            parts = ID.split("_from_")
            num = parts[0]
            setattr(DEN, parts[1], cms.vstring("pass"))
        if scenario.find("tagiso") != -1:  
            DEN.tag_combRelIso = cms.vdouble(-1, 0.1)
        if scenario.find("probeiso") != -1:
            DEN.isoTrk03Abs = cms.vdouble(-1, 3)
        #if scenario.find("calo") != -1: DEN.caloCompatibility = cms.vdouble(0.9,1.1)  # same as above, I think.
        setattr(module.Efficiencies, ID+"_"+X, cms.PSet(
            EfficiencyCategoryAndState = cms.vstring(num,"pass"),
            UnbinnedVariables = cms.vstring("mass"),
            BinnedVariables = DEN,
            BinToPDFmap = cms.vstring(shape)
        ))
        if scenario.find("mc") != -1:
            setattr(module.Efficiencies, ID+"_"+X+"_mcTrue", cms.PSet(
                EfficiencyCategoryAndState = cms.vstring(num,"pass"),
                UnbinnedVariables = cms.vstring("mass"),
                BinnedVariables = DEN.clone(mcTrue = cms.vstring("true"))
            ))
        setattr(process, "TnP_MuonID_"+ID+"_"+X, module)        
        setattr(process, "run_"+ID+"_"+X, cms.Path(module))


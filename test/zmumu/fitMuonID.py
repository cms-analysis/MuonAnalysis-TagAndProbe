import FWCore.ParameterSet.Config as cms
### USAGE:
###    cmsRun fitMuonID.py <scenario> [ <id> [ <binning1> ... <binningN> ] ]
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

    Variables = cms.PSet(
        mass = cms.vstring("Tag-muon Mass", "76", "125", "GeV/c^{2}"),
        pt = cms.vstring("muon p_{T}", "0", "1000", "GeV/c"),
        eta    = cms.vstring("muon #eta", "-2.5", "2.5", ""),
        abseta = cms.vstring("muon |#eta|", "0", "2.5", ""),
        tag_nVertices = cms.vstring("Number of vertices", "0", "999", ""),
        pfCombRelIso04EACorr = cms.vstring("Number of vertices", "0", "999", ""),
        SIP = cms.vstring("Number of vertices", "0", "999", ""),
        run = cms.vstring("Number of vertices", "-999", "999999", ""),
        HLT_Mu8  = cms.vstring("Mu8 leg",  "-1", "2", ""),
        HLT_Mu17 = cms.vstring("Mu17 leg", "-1", "2", ""),
    ),

    Categories = cms.PSet(
        Glb   = cms.vstring("Global", "dummy[pass=1,fail=0]"),
        PF    = cms.vstring("PF Muon", "dummy[pass=1,fail=0]"),
        TM    = cms.vstring("Tracker Muon", "dummy[pass=1,fail=0]"),
        mvaIsoCut = cms.vstring("MC true", "dummy[pass=1,fail=0]"),
        mcTrue = cms.vstring("MC true", "dummy[true=1,false=0]"),
    ),
    Cuts = cms.PSet(),
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
        ),
        vpvPlusExpoMin70 = cms.vstring(
            "Voigtian::signal1(mass, mean1[90,80,100], width[2.495], sigma1[2,1,3])",
            "Voigtian::signal2(mass, mean2[90,80,100], width,        sigma2[4,3,10])",
            "SUM::signal(vFrac[0.8,0.5,1]*signal1, signal2)",
            "Exponential::backgroundPass(mass, lp[-0.1,-1,0.1])",
            "Exponential::backgroundFail(mass, lf[-0.1,-1,0.1])",
            "efficiency[0.9,0.7,1]",
            "signalFractionInPassing[0.9]"
        )
    ),

    binnedFit = cms.bool(True),
    binsForFit = cms.uint32(40),
    saveDistributionsPlot = cms.bool(False),

    Efficiencies = cms.PSet(), # will be filled later
)

TRIGGER = cms.PSet(tag_Mu24 = cms.vstring("pass"))
if "mc" in scenario or "39X" in scenario or "38X" in scenario:
    TRIGGER = cms.PSet(tag_Mu15 = cms.vstring("pass"), tag_pt = cms.vdouble(24.,9999.))

# PT_ETA_BINS = cms.PSet(
#     pt     = cms.vdouble(  3, 5, 7, 10, 15, 20, 30, 40, 100 ),
#     abseta = cms.vdouble(  0.0, 1.2, 2.4)
# )

# Cris' bins
PT_ETA_BINS = cms.PSet(
    pt     = cms.vdouble(  10, 20, 30, 40, 60, 100 ),
    abseta = cms.vdouble(  0.0, 1.2, 2.4)
)


PT_ETA_BINS_MU17 = cms.PSet(
    pt     = cms.vdouble(  10, 15, 16, 16.5, 17, 17.5, 18.0, 19, 20, 30, 40, 100 ),
    abseta = cms.vdouble(  0.0, 1.2, 2.4)
)
PT_ETA_BINS_MU8  = cms.PSet(
    pt     = cms.vdouble(  3, 5, 6, 7, 7.5, 8, 8.5, 9, 10, 12, 15, 20, 30, 40, 100 ),
    abseta = cms.vdouble(  0.0, 1.2, 2.4)
)

# ETA_BINS = cms.PSet(
#     pt  = cms.vdouble(20,100),
#     eta = cms.vdouble(*[x/10. for x in xrange(-24,25,1)]),    
# )


# Cris' bins
ETA_BINS = cms.PSet(
    pt  = cms.vdouble(20,100),
    eta = cms.vdouble(-2.4, -2.1, -1.6, -1.2, -0.9, -0.4, 0., 0.4, 0.9, 1.2, 1.6, 2.1, 2.4),
    #eta = cms.vdouble(-2.4, -2.1, -1.6, -1.1, -0.6, 0, 0.6, 1.1, 1.6, 2.1, 2.4),
)



# VTX_BINS  = cms.PSet(
#     pt     = cms.vdouble(  20, 120 ),
#     abseta = cms.vdouble(  0.0, 2.4),
#     tag_nVertices = cms.vdouble(*[x-0.5 for x in range(21)])
# )


# Cris' bins
VTX_BINS  = cms.PSet(
    pt     = cms.vdouble(  20, 120 ),
    abseta = cms.vdouble(  0.0, 2.4),
    #tag_nVertices = cms.vdouble(0.5,2.5,4.5,6.5,8.5,10.5,12.5,14.5,16.5,18.5,20.5)
    tag_nVertices = cms.vdouble(0.5,2.5,4.5,6.5,8.5,10.5,12.5,16.5,20.5,25,30,40)
)


PREFIX="/afs/cern.ch/work/g/gpetrucc/CMSSW_5_2_4_patch4/src/MuonAnalysis/TagAndProbe/test/zmumu/"

process.TnP_MuonID = Template.clone(
    InputFileNames = cms.vstring(),
    InputTreeName = cms.string("fitter_tree"),
    InputDirectoryName = cms.string("tpTree"),
    OutputFileName = cms.string("TnP_MuonID_%s.root" % scenario),
    Efficiencies = cms.PSet(),
)
if "_weight" in scenario:
    process.TnP_MuonID.WeightVariable = cms.string("weight")
    process.TnP_MuonID.Variables.weight = cms.vstring("weight","0","10","")


if "data2011" in scenario:
    process.TnP_MuonID.InputFileNames = [ PREFIX+"tnpZ_Run2011_forID.root" ]
if "mc2011" in scenario:
    process.TnP_MuonID.InputFileNames = [ PREFIX+"tnpZ_MC_DY50Fall11_forID_withNVtxWeights2011.root" ]
if "data2012" in scenario:
    process.TnP_MuonID.InputFileNames = [ PREFIX+"tnpZ_Data.190456-193557_forID.root" ]
if "mc2012" in scenario:
    process.TnP_MuonID.InputFileNames = [ PREFIX+"tnpZ_MC_DYJetsToLL_S7_START52_V9_forID_withNVtxWeights190456-193557.root" ]

#IDS = [ "PF", "SIP4_from_PF", "PFIso25_from_PF_and_SIP4", "PFIso40_from_PF_and_SIP4", "mvaIsoCut_from_PF_and_SIP4" ]
IDS = [ "PF", "SIP4_from_PF", "PFIso40_from_PF_and_SIP4"]
ALLBINS = [("pt_abseta",PT_ETA_BINS),("eta", ETA_BINS)]
ALLBINS += [ ("vtx",VTX_BINS)]

if len(args) > 1 and args[1] not in IDS: IDS += [ args[1] ]
for ID in IDS:
    if len(args) > 1 and ID != args[1]: continue
    for X,B in ALLBINS:
        if len(args) > 2 and X not in args[2:]: continue
        module = process.TnP_MuonID.clone(OutputFileName = cms.string("TnP_MuonID_%s_%s_%s.root" % (scenario, ID, X)))
        shape = "vpvPlusExpo"
        #if "eta" in X and not "abseta" in X: shape = "voigtPlusExpo"
        #if "pt_abseta" in X: shape = "voigtPlusExpo"
        if "Mu17" in ID and "pt_abseta" in X: B = PT_ETA_BINS_MU17
        if "Mu8" in ID and "pt_abseta" in X: B = PT_ETA_BINS_MU8
        DEN = B.clone(); num = ID; numstate = "pass"
        if "_from_" in ID:
            parts = ID.split("_from_")
            num = parts[0]
            for D in parts[1].split("_and_"):
                if D == "SIP4":      DEN.SIP = cms.vdouble(0,4.0)
                elif D == "PFIso25": DEN.pfCombRelIso04EACorr  = cms.vdouble(0,0.25)
                elif D == "PFIso40": DEN.pfCombRelIso04EACorr  = cms.vdouble(0,0.40)
                else: setattr(DEN, D, cms.vstring("pass"))
        numString = cms.vstring()
        for N in num.split("_and_"):
            if N == "SIP4":
                module.Cuts.SIP4 = cms.vstring("SIP4","SIP", "4.0") 
                numString += [N, "below"]
            elif N == "PFIso25":
                module.Cuts.PFIso25 = cms.vstring("PFIso","pfCombRelIso04EACorr", "0.25") 
                numString += [N, "below"]
            elif N == "PFIso40":
                module.Cuts.PFIso40 = cms.vstring("PFIso","pfCombRelIso04EACorr", "0.40") 
                numString += [N, "below"]
            elif N in [ "HLT_Mu17", "HLT_Mu8" ]:
                setattr(module.Cuts, N+"_Pass", cms.vstring(N+"_Pass", N, "0.5")) 
                numString += [N+"_Pass", "above"]
            else:
                numString += [N, "pass"]
        setattr(module.Efficiencies, ID+"_"+X, cms.PSet(
            EfficiencyCategoryAndState = numString,
            UnbinnedVariables = cms.vstring("mass"),
            BinnedVariables = DEN,
            BinToPDFmap = cms.vstring(shape)
        ))
        if scenario.find("mc") != -1:
            setattr(module.Efficiencies, ID+"_"+X+"_mcTrue", cms.PSet(
                EfficiencyCategoryAndState = numString,
                UnbinnedVariables = cms.vstring("mass"),
                BinnedVariables = DEN.clone(mcTrue = cms.vstring("true"))
            ))
            if "_weight" in scenario:
                getattr(module.Efficiencies, ID+"_"+X          ).UnbinnedVariables.append("weight")
                getattr(module.Efficiencies, ID+"_"+X+"_mcTrue").UnbinnedVariables.append("weight")
        setattr(process, "TnP_MuonID_"+ID+"_"+X, module)        
        setattr(process, "run_"+ID+"_"+X, cms.Path(module))


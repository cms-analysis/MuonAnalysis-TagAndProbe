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
        tkTrackerLay = cms.vstring("number of valid tracker layers", "-0.5", "30.5", ""),
        tkPixelLay   = cms.vstring("number of valid pixel layers", "-0.5", "12.5", ""),
        numberOfMatchedStations = cms.vstring("number of matched stations", "-0.5", "7.5", ""),
        tag_nVertices = cms.vstring("Number of vertices", "0", "999", ""),
        pfCombRelIso04EACorr = cms.vstring("Number of vertices", "0", "999", ""),
        combRelIsoPF04dBeta  = cms.vstring("Number of vertices", "0", "999", ""),
        SIP = cms.vstring("Number of vertices", "0", "999", ""),
        run = cms.vstring("Number of vertices", "-999", "999999", ""),
        dzPV = cms.vstring("Number of vertices", "-999", "999", ""),
        dB   = cms.vstring("Number of vertices", "-999", "999", ""),
    ),

    Categories = cms.PSet(
        Glb   = cms.vstring("Global", "dummy[pass=1,fail=0]"),
        GlbPT = cms.vstring("Global", "dummy[pass=1,fail=0]"),
        PF    = cms.vstring("PF Muon", "dummy[pass=1,fail=0]"),
        TMA   = cms.vstring("Tracker Muon", "dummy[pass=1,fail=0]"),
        Tight2012 = cms.vstring("Tight",   "dummy[pass=1,fail=0]"),
        mcTrue = cms.vstring("MC true", "dummy[true=1,false=0]"),
        HLT_Mu17 = cms.vstring("Mu17", "dummy[pass=1,fail=0]"),
        HLT_Mu8 = cms.vstring("Mu8", "dummy[pass=1,fail=0]"),
        HLT_TkMu17 = cms.vstring("TkMu17", "dummy[pass=1,fail=0]"),
        HLT_TkMu8 = cms.vstring("TkMu8", "dummy[pass=1,fail=0]"),
        HLT_OrMu17 = cms.vstring("OrMu17", "dummy[pass=1,fail=0]"),
        HLT_OrMu8 = cms.vstring("OrMu8", "dummy[pass=1,fail=0]"),
        DoubleMu17Mu8_Mu17 = cms.vstring("Mu17", "dummy[pass=1,fail=0]"),
        DoubleMu17Mu8_Mu8 = cms.vstring("Mu8", "dummy[pass=1,fail=0]"),
        DoubleMu17TkMu8_Mu17 = cms.vstring("TkMu17", "dummy[pass=1,fail=0]"),
        DoubleMu17TkMu8_TkMu8 = cms.vstring("TkMu8", "dummy[pass=1,fail=0]"),
    ),
    Expressions = cms.PSet(
       ZZLooseVar = cms.vstring("", "Glb==1 || numberOfMatchedStations > 0", "Glb", "numberOfMatchedStations"),
       Tight2012noIPVar = cms.vstring("", "GlbPT==1 && numberOfMatchedStations > 1 && tkTrackerLay > 5 && tkPixelLay > 0", 
                    "GlbPT", "numberOfMatchedStations", "tkTrackerLay", "tkPixelLay"),
       Tight2012withIPVar = cms.vstring("", "Tight2012==1 && abs(dzPV) < 0.5", 
                    "Tight2012", "dzPV"),
    ),
    Cuts = cms.PSet(
      ZZLoose         = cms.vstring("ZZLoose",  "ZZLooseVar",   "0.5"),
      Tight2012noIP   = cms.vstring("Tight2012noIP",  "Tight2012noIPVar",   "0.5"),
      Tight2012withIP = cms.vstring("Tight2012withIP","Tight2012withIPVar", "0.5"),
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

PT_ETA_BINS = cms.PSet(
    pt     = cms.vdouble(  3, 5, 7, 10, 15, 20, 30, 40, 100 ),
    abseta = cms.vdouble(  0.0, 1.2, 2.4)
)

PT_ETA_BINS2 = cms.PSet(
    pt     = cms.vdouble(  3, 5, 10, 15, 20, 30, 40, 100 ),
    abseta = cms.vdouble(  0.0, 1.2, 2.4 )
)


PT_ETA_BINS_MU17 = cms.PSet(
    pt     = cms.vdouble(  10, 15, 16, 12.25, 16.5, 16.75, 17, 17.25, 17.5, 18.75, 18.0, 19, 20, 30, 40, 100 ),
    abseta = cms.vdouble(  0.0, 1.2, 2.4)
)
PT_ETA_BINS_MU8  = cms.PSet(
    pt     = cms.vdouble(  3, 5, 6, 7, 7.5, 7.75, 8, 8.25, 8.5, 8.75, 9, 10, 12, 15, 20, 30, 40, 100 ),
    abseta = cms.vdouble(  0.0, 1.2, 2.4)
)

# Ivan's bins
ETA_BINS = cms.PSet(
    pt  = cms.vdouble(20,100),
    eta = cms.vdouble(-2.4, -2.1, -1.6, -1.2, -0.9, -0.6, -0.3, -0.2, 0.2, 0.3, 0.6, 0.9, 1.2, 1.6, 2.1, 2.4),
    #eta = cms.vdouble(-2.4, -2.1, -1.6, -0.9, -0.4, 0, 0.4, 0.9, 1.6, 2.1, 2.4),
)


# Cris' bins
VTX_BINS  = cms.PSet(
    pt     = cms.vdouble(  20, 120 ),
    abseta = cms.vdouble(  0.0, 2.4),
    #tag_nVertices = cms.vdouble(0.5,2.5,4.5,6.5,8.5,10.5,12.5,14.5,16.5,18.5,20.5)
    tag_nVertices = cms.vdouble(0.5,2.5,4.5,6.5,8.5,10.5,12.5,16.5,20.5,25,30,40)
    #tag_nVertices = cms.vdouble(0.5,4.5,6.5,8.5,10.5,12.5,16.5,20.5,25,30)
)


PREFIX="/afs/cern.ch/work/g/gpetrucc/CMSSW_5_2_4_patch4/src/MuonAnalysis/TagAndProbe/test/zmumu/ZZ4L/Cris/"
if "2012" in scenario:
    PREFIX="/afs/cern.ch/work/g/gpetrucc/CMSSW_5_2_4_patch4/src/MuonAnalysis/TagAndProbe/test/zmumu/ZZ4L/Cris_5fb/"

process.TnP_MuonID = Template.clone(
    InputFileNames = cms.vstring(),
    InputTreeName = cms.string("fitter_tree"),
    InputDirectoryName = cms.string("tpTree"),
    OutputFileName = cms.string("TnPZ_MuonID_%s.root" % scenario),
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
    process.TnP_MuonID.InputFileNames = [ PREFIX+"tnpZ_Data_Run2012_forID.root" ]
if "data2012C" in scenario:
    process.TnP_MuonID.InputFileNames = [ PREFIX+"tnpZ_Data_Run2012C_upTo199011_forID.root" ]
if "mc2012" in scenario:
    process.TnP_MuonID.InputFileNames = [ PREFIX+"tnpZ_MC_DYJetsToLL_madgraph_Summer12_PU_S7_START52_V9_v2_forID_withNVtxWeights190456-196531.root" ]
if "mc2012C" in scenario:
    process.TnP_MuonID.InputFileNames = [ PREFIX+"tnpZ_MC_DYJetsToLL_forID_withNVtxWeights_Run2012C_upTo199011.root" ]

IDS = [ "PF", "SIP4_from_PF", "PFIso40_from_PF_and_SIP4" ]
ALLBINS = [("pt_abseta",PT_ETA_BINS),("eta", ETA_BINS)]
ALLBINS += [ ("vtx",VTX_BINS)]
ALLBINS += [("pt_abseta2",PT_ETA_BINS2) ]

if len(args) > 1 and args[1] not in IDS: IDS += [ args[1] ]
for ID in IDS:
    if len(args) > 1 and ID != args[1]: continue
    for X,B in ALLBINS:
        if len(args) > 2 and X not in args[2:]: continue
        module = process.TnP_MuonID.clone(OutputFileName = cms.string("TnP_MuonID_%s_%s_%s.root" % (scenario, ID, X)))
        shape = "vpvPlusExpo"
        #if "eta" in X and not "abseta" in X: shape = "voigtPlusExpo"
        #if "pt_abseta" in X: shape = "voigtPlusExpo"
        if "Mu8" in ID and "pt_abseta" in X: B = PT_ETA_BINS_MU8
        if ("_Mu17" in ID or "_OrMu17" in ID or "_TkMu17" in ID) and "pt_abseta" in X: B = PT_ETA_BINS_MU17 ## DONT's swap with above
        DEN = B.clone(); num = ID; numstate = "pass"
        if "_from_" in ID:
            parts = ID.split("_from_")
            num = parts[0]
            for D in parts[1].split("_and_"):
                if D == "SIP4only": 
                    DEN.SIP = cms.vdouble(0,4.0)
                elif D == "SIP4": 
                    DEN.SIP = cms.vdouble(0,4.0)
                    DEN.dB = cms.vdouble(0,0.5)
                    DEN.dzPV = cms.vdouble(-1.0,1.0)
                elif D == "Tight2012noIP":
                    DEN.GlbPT = cms.vstring("pass")
                    DEN.numberOfMatchedStations = cms.vdouble(1.5,99)
                    DEN.tkTrackerLay = cms.vdouble(5.5,99)
                    DEN.tkPixelLay = cms.vdouble(0.5,99)
                elif D == "Tight2012withIP":
                    DEN.Tight2012 = cms.vstring("pass")
                    DEN.dzPV = cms.vdouble(-0.5,0.5)
                elif D == "PFIso40": DEN.pfCombRelIso04EACorr  = cms.vdouble(0,0.40)
                elif D == "PFIso12DB": DEN.combRelIsoPF04dBeta = cms.vdouble(0,0.12)
                elif D == "PFIso20DB": DEN.combRelIsoPF04dBeta = cms.vdouble(0,0.20)
                elif D == "PostTS" : DEN.run = cms.vdouble(192000,999999)
                elif D == "2011A" : DEN.run = cms.vdouble(160404,173692)
                elif D == "2011B" : DEN.run = cms.vdouble(175832,177791)
                else: setattr(DEN, D, cms.vstring("pass"))
        numString = cms.vstring()
        for N in num.split("_and_"):
            if N == "SIP4":
                module.Cuts.SIP4 = cms.vstring("SIP4","SIP", "4.0") 
                module.Cuts.dB05 = cms.vstring("dB05","dB",  "0.5") 
                module.Cuts.dZp1 = cms.vstring("dZp1","dzPV",  "1.0") 
                module.Cuts.dZm1 = cms.vstring("dZm1","dzPV", "-1.0") 
                numString += [N, "below", "dB05", "below", "dZp1", "below", "dZm1", "above" ]
            elif N == "SIP4only":
                module.Cuts.SIP4 = cms.vstring("SIP4","SIP", "4.0") 
                numString += ["SIP4", "below"]
            elif N == "PFIso40":
                module.Cuts.PFIso40 = cms.vstring("PFIso","pfCombRelIso04EACorr", "0.40") 
                numString += [N, "below"]
            elif N == "PFIso20DB":
                module.Cuts.PFIso20DB = cms.vstring("PFIso","combRelIsoPF04dBeta", "0.20") 
                numString += [N, "below"]
            elif N == "PFIso12DB":
                module.Cuts.PFIso12DB = cms.vstring("PFIso","combRelIsoPF04dBeta", "0.12") 
                numString += [N, "below"]
            elif N in [ "Tight2012noIP",  "Tight2012withIP", "ZZLoose" ]:
                numString += [N, "above" ]
            else:
                numString += [N, "pass" ]
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


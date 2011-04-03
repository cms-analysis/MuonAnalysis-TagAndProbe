import FWCore.ParameterSet.Config as cms

### USAGE:
###    cmsRun fitMuonTrigger.py <scenario>
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

process.TnP_MuonTrigger = cms.EDAnalyzer("TagProbeFitTreeAnalyzer",
                                         NumCPU = cms.uint32(1),
                                         SaveWorkspace = cms.bool(True),

                                         InputFileNames = cms.vstring('rfio:/castor/cern.ch/user/t/trocino/Temp/tnpJPsi_Data_Nov4B_L3.root'
                                                                      ),
                                         InputTreeName = cms.string("fitter_tree"),
                                         InputDirectoryName = cms.string("tpTree"),
                                         OutputFileName = cms.string("L3Efficiency_%s_Jpsi.root" % scenario),
                                         #                                     InputFileNames = cms.vstring('rfio:/castor/cern.ch/user/b/botta/TnP_Paper_Tree/tnpJPsi_MuOnia_Run2010B_Nov4.root'),
                                         #                                     InputTreeName = cms.string("fitter_tree"),
                                         #                                     InputDirectoryName = cms.string("tpTree"),
                                         #                                     OutputFileName = cms.string("/tmp/botta/L3Efficiency_%s_Jpsi.root" % scenario),

                                         Variables = cms.PSet(mass   = cms.vstring("Tag-Probe Mass", "2.8", "3.4", "GeV/c^{2}"),
                                                              pt     = cms.vstring("Probe p_{T}", "0", "1000", "GeV/c"),
                                                              eta    = cms.vstring("Probe #eta", "-2.5", "2.5", ""),
                                                              abseta = cms.vstring("Probe |#eta|", "0.", "2.5", ""),
                                                              phi    = cms.vstring("Probe #phi", "-3.15.", "3.15", ""),
                                                              tag_pt = cms.vstring("Tag p_{T}", "2.6", "1000", "GeV/c"),
                                                              run = cms.vstring("RunNumber","133446","149442", "RunNumber"),
                                                              pair_dphiVtxTimesQ = cms.vstring("pair_dphiVtxTimesQ", "-1", "1", "deltaPhi*Charge"),
                                                              pair_distM1  = cms.vstring("q1 * (#phi1-#phi2)", "-99999", "999999", ""),
                                                              #charge = ##################
                                                              l1q    = cms.vstring("L1 matched object quality", "-1000", "10", ""),
                                                              l1pt   = cms.vstring("L1 matched object p_{T}", "0", "1000", "GeV/c"),
                                                              l2pt   = cms.vstring("L2 matched object p_{T}", "0", "1000", "GeV/c"),
                                                              l3pt   = cms.vstring("L3 matched object p_{T}", "0", "1000", "GeV/c"),
                                                         
                                                              ),


                                         Categories = cms.PSet(
        VBTF      = cms.vstring("VBTF-like #mu",   "dummy[pass=1,fail=0]"), 
        Mu9       = cms.vstring("HLTMu9",             "dummy[pass=1,fail=0]"),
        Mu15      = cms.vstring("HLTMu15",            "dummy[pass=1,fail=0]"),
        
        TMOST     = cms.vstring("TMOST",          "dummy[pass=1,fail=0]"),
        
        tag_Mu0_L2Mu0_MU = cms.vstring("TagTrigger_Mu0", "dummy[pass=1,fail=0]"),
        tag_Mu3_L2Mu0_MU = cms.vstring("TagTrigger_Mu3", "dummy[pass=1,fail=0]"),
        tag_Mu5_L2Mu0_MU = cms.vstring("TagTrigger_Mu5", "dummy[pass=1,fail=0]"),
        MuX_L2Mu0_L2 = cms.vstring("ProbeTrigger", "dummy[pass=1,fail=0]"),
        
        DoubleMu0_Quarkonium = cms.vstring("PassingProbeTrigger", "dummy[pass=1,fail=0]"),
        DoubleMu0 = cms.vstring("PassingProbeTrigger", "dummy[pass=1,fail=0]"),
        
        mcTrue = cms.vstring("MC true", "dummy[true=1,false=0]"),
        ),

                                         
                                         Expressions = cms.PSet(OrMu9Mu15 = cms.vstring("OrMu9Mu15", "Mu15 == 1 || ( Mu9 == 1 && l3pt >= 15 )", "Mu15", "Mu9", "l3pt"),
                                                                OrDoubleMu0 = cms.vstring("OrDoubleMu0", "DoubleMu0 ==1 || DoubleMu0_Quarkonium==1", "DoubleMu0", "DoubleMu0_Quarkonium"),
                                                                OrDoubleMu3 = cms.vstring("OrDoubleMu0", "( DoubleMu0 ==1 && l3pt >= 3) || ( DoubleMu0_Quarkonium==1 && l3pt >= 3)", "DoubleMu0", "DoubleMu0_Quarkonium", "l3pt"),
                                                                
                                                                OrDoubleMu9 = cms.vstring("OrDoubleMu0", "( DoubleMu0 ==1 && l3pt >= 9) || ( DoubleMu0_Quarkonium==1 && l3pt >= 9)", "DoubleMu0", "DoubleMu0_Quarkonium", "l3pt"),
                                                                #OrMuXL2Mu0MU = cms.vstring("OrMuXL2Mu0MU", "(tag_Mu0_L2Mu0_MU ==1 || tag_Mu3_L2Mu0_MU==1 || tag_Mu5_L2Mu0_MU==1) && MuX_L2Mu0_L2==1 ", "tag_Mu0_L2Mu0_MU", "tag_Mu3_L2Mu0_MU", "tag_Mu5_L2Mu0_MU", "MuX_L2Mu0_L2==1"),
                                                                OrDoubleMu12 = cms.vstring("OrDoubleMu0", "( DoubleMu0 ==1 && l3pt >= 12) || ( DoubleMu0_Quarkonium==1 && l3pt >= 12)", "DoubleMu0", "DoubleMu0_Quarkonium", "l3pt"),
                                                                ),
                                         
                                         Cuts = cms.PSet( L1Mu0Pt  = cms.vstring("L1Mu3Pt",  "l1pt", "0"),
                                                          L1Mu7Pt  = cms.vstring("L1Mu7Pt",   "l1pt", "7.01"),
                                                          L1Mu0Q   = cms.vstring("L1Mu0Q",   "l1q",  "1.5"),
                                                          L1Mu7Q   = cms.vstring("L1Mu7Q",  "l1q",  "3.5"),
                                                          L2Mu7Pt  = cms.vstring("L1Mu15Pt", "l2pt", "7.01"),
                                                          #RunWithMu15  = cms.vstring("RunWithMu15", "run",  "147180"),
                                                          
                                                          SelOrMu9Mu15 = cms.vstring("SelOrMu9Mu15", "OrMu9Mu15", "0.5"),
                                                          SelOrDoubleMu0 = cms.vstring("SelOrDoubleMu0", "OrDoubleMu0", "0.5"),
                                                          SelOrDoubleMu3 = cms.vstring("SelOrDoubleMu3", "OrDoubleMu3", "0.5"),
                                                          SelOrDoubleMu9 = cms.vstring("SelOrDoubleMu9", "OrDoubleMu9", "0.5"),
                                                          SelOrDoubleMu12 = cms.vstring("SelOrDoubleMu12", "OrDoubleMu12", "0.5"),
                                                          #SelOrMuXL2Mu0MU = cms.vstring("SelOrMuXL2Mu0MU", "OrMuXL2Mu0MU", "0.5"),
                                                          ),
                                         
                                         PDFs = cms.PSet(gaussPlusExpo = cms.vstring("CBShape::signal(mass, mean[3.1,3.0,3.2], sigma[0.05,0.02,0.1], alpha[3., 0., 4.], n[1, 0., 100.])",
                                                                                     "Chebychev::backgroundPass(mass, {cPass[0,-1,1], cPass2[0,-1,1]})",
                                                                                     "Chebychev::backgroundFail(mass, {cFail[0,-1,1], cFail2[0,-1,1]})",
                                                                                     "efficiency[0.9,0,1]",
                                                                                     "signalFractionInPassing[0.9]"
                                                                                     )
                                                         ),
                                         
                                         binnedFit = cms.bool(True),
                                         binsForFit = cms.uint32(40),
                                         Efficiencies = cms.PSet(), # will be filled later
                                         )


#if scenario == "datalike_mc":
#   process.TnP_MuonTrigger.InputFileNames = [ "tnpZ_MC.root", ]

#################
# HLT_Mu0 vs L2 #   
# HLT_Mu3 vs L2 #  
#################
## BIN DEFINITIONS FOR THE J/PSI T&P EFFICIENCIES MEASUREMENTS

ALLBINS=[ ("pt_eta",
           cms.PSet(pt     = cms.vdouble(0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.75, 9.25, 9.75, 10.25, 11.0, 14.0, 17.0, 20.0),
                    abseta = cms.vdouble(0.0, 0.9, 2.1, 2.4),
                    pair_dphiVtxTimesQ = cms.vdouble(-2, 0),
                    pair_distM1 = cms.vdouble(200,1000),
                    )
           ),
          ("pt",
           cms.PSet(pt     = cms.vdouble(0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.75, 9.25, 9.75, 10.25, 11.0, 14.0, 17.0, 20.0),
                    abseta = cms.vdouble(0.0, 2.4),
                    pair_dphiVtxTimesQ = cms.vdouble(-2, 0),
                    pair_distM1 = cms.vdouble(200,1000),
                    )
           ),
          ("eta",
           cms.PSet(eta  = cms.vdouble(-2.4, -2.1, -1.6, -1.2, -0.9, -0.4, 0., 0.4, 0.9, 1.2, 1.6, 2.1, 2.4),
                    pair_dphiVtxTimesQ = cms.vdouble(-2, 0),
                    pair_distM1 = cms.vdouble(200,1000),
                    pt = cms.vdouble(3., 20.0),                    
                    )
           )
          ]

for (T,Tn,M) in [ ("SelOrDoubleMu9","L3Mu9overL2Mu7","VBTF")]:
    for BN,BV in ALLBINS:
        BINNEDVARS = BV.clone()

        ### already asked with the subTree
        #setattr(BINNEDVARS, "SelOrMuXL2Mu0MU", cms.vstring("above"))
        #setattr(BINNEDVARS, "MuX_L2Mu0_L2", cms.vstring("pass"))
        setattr(BINNEDVARS, M, cms.vstring("pass"))


        if T == "SelOrDoubleMu9":
            if BN == "eta":
                BINNEDVARSMOD = BINNEDVARS.clone(pt = cms.vdouble(11.0, 20.0),
                                                 l1q  = cms.vdouble(3.5, 10.),
                                                 l1pt = cms.vdouble(7.01, 1000.),
                                                 l2pt = cms.vdouble(7.01, 1000.),
                                                 )
            else:
                BINNEDVARSMOD = BINNEDVARS.clone(l1q  = cms.vdouble(3.5, 10.),
                                                 l1pt = cms.vdouble(7.01, 1000.),
                                                 l2pt = cms.vdouble(7.01, 1000.),
                                                 )
            setattr(process.TnP_MuonTrigger.Efficiencies, M+"_To_"+Tn+"_"+BN, cms.PSet(
                EfficiencyCategoryAndState = cms.vstring(T,"above"), # above because T is a Cut
                UnbinnedVariables = cms.vstring("mass"),
                BinnedVariables = BINNEDVARSMOD,
                BinToPDFmap = cms.vstring("gaussPlusExpo")
                ))


process.p = cms.Path(process.TnP_MuonTrigger)

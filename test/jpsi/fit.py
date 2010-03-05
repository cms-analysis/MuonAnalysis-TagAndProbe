PTBINS = cms.vdouble(3.5, 4.5, 50.0)
ETABINS = cms.vdouble(-2.1, 0.0, 2.1)
FILEPREFIX = "coarse_"

import FWCore.ParameterSet.Config as cms

process = cms.Process("TagProbe")

process.load('FWCore.MessageService.MessageLogger_cfi')

process.source = cms.Source("EmptySource")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

Template = cms.EDAnalyzer("TagProbeFitTreeAnalyzer",
    NumCPU = cms.uint32(8),
    SaveWorkspace = cms.bool(False),

    Variables = cms.PSet(
        mass = cms.vstring("Tag-Probe Mass", "2.8", "3.5", "GeV/c^{2}"),
        pt = cms.vstring("Probe p_{T}", "0", "1000", "GeV/c"),
        p = cms.vstring("Probe p", "0", "1000", "GeV/c"),
        eta = cms.vstring("Probe #eta", "-2.5", "2.5", ""),
        phi = cms.vstring("Probe #phi", "-3.1416", "3.1416", "")
    ),

    Categories = cms.PSet(
        mcTrue = cms.vstring("MC true", "dummy[true=1,false=0]"),
        passingGlb = cms.vstring("isGlobalMuon", "dummy[pass=1,fail=0]"),
        passingTrkEx = cms.vstring("isTrackerMuon", "dummy[pass=1,fail=0]"),
        Glb = cms.vstring("isGlobalMuon", "dummy[pass=1,fail=0]"),
        Trk = cms.vstring("isTrackerMuon", "dummy[pass=1,fail=0]"),
        HLTMu3 = cms.vstring("isHLT_Mu3", "dummy[pass=1,fail=0]"),
        L1DiMuOpen = cms.vstring("isHLT_L1DoubleMuOpen", "dummy[pass=1,fail=0]")
    ),

    PDFs = cms.PSet(
        gaussPlusLinear = cms.vstring(
            "Gaussian::signal(mass, mean[3.1,3.0,3.2], sigma[0.03,0.01,0.05])",
            "Chebychev::backgroundPass(mass, cPass[0,-1,1])",
            "Chebychev::backgroundFail(mass, cFail[0,-1,1])",
            "efficiency[0.9,0,1]",
            "signalFractionInPassing[0.9]"
        ),
        gaussPlusQuadratic = cms.vstring(
            "Gaussian::signal(mass, mean[3.1,3.0,3.2], sigma[0.03,0.01,0.05])",
            "Chebychev::backgroundPass(mass, {cPass1[0,-1,1], cPass2[0,-1,1]})",
            "Chebychev::backgroundFail(mass, {cFail1[0,-1,1], cFail2[0,-1,1]})",
            "efficiency[0.9,0,1]",
            "signalFractionInPassing[0.9]"
        )
    )
)

process.TnP_MuFromTk = Template.clone(
    InputFileNames = cms.vstring("tnpJPsi_MuonID_JPsi-0.084pb.root", "tnpJPsi_MuonID_ppMuX.root"),
    InputDirectoryName = cms.string("histoMuFromTk"),
    InputTreeName = cms.string("fitter_tree"),
    OutputFileName = cms.string(FILEPREFIX+"TnP_MuFromTk.root"),
    Efficiencies = cms.PSet(
        Glb_pt_eta = cms.PSet(
            EfficiencyCategoryAndState = cms.vstring("Glb","pass"),
            UnbinnedVariables = cms.vstring("mass"),
            BinnedVariables = cms.PSet(
                pt = PTBINS,
                eta = ETABINS
            ),
            BinToPDFmap = cms.vstring("gaussPlusLinear")
        ),
        Glb_pt_eta_mcTrue = cms.PSet(
            EfficiencyCategoryAndState = cms.vstring("Glb","pass"),
            UnbinnedVariables = cms.vstring("mass"),
            BinnedVariables = cms.PSet(
                mcTrue = cms.vstring("true"),
                pt = PTBINS,
                eta = ETABINS
            ),
            BinToPDFmap = cms.vstring()
        ),
        Trk_pt_eta = cms.PSet(
            EfficiencyCategoryAndState = cms.vstring("Trk","pass"),
            UnbinnedVariables = cms.vstring("mass"),
            BinnedVariables = cms.PSet(
                pt = PTBINS,
                eta = ETABINS
            ),
            BinToPDFmap = cms.vstring("gaussPlusLinear")
        ),
        Trk_pt_eta_mcTrue = cms.PSet(
            EfficiencyCategoryAndState = cms.vstring("Trk","pass"),
            UnbinnedVariables = cms.vstring("mass"),
            BinnedVariables = cms.PSet(
                mcTrue = cms.vstring("true"),
                pt = PTBINS,
                eta = ETABINS
            ),
            BinToPDFmap = cms.vstring()
        ),
    )
)

process.TnP_MuFromCal = Template.clone(
    InputFileNames = cms.vstring("tnpJPsi_MuonID_JPsi-0.084pb.root", "tnpJPsi_MuonID_ppMuX.root"),
    InputDirectoryName = cms.string("histoMuFromCal"),
    InputTreeName = cms.string("fitter_tree"),
    OutputFileName = cms.string(FILEPREFIX+"TnP_MuFromCal.root"),
    Efficiencies = cms.PSet(
        Glb_pt_eta = cms.PSet(
            EfficiencyCategoryAndState = cms.vstring("passingGlb","pass"),
            UnbinnedVariables = cms.vstring("mass"),
            BinnedVariables = cms.PSet(
                pt = PTBINS,
                eta = ETABINS
            ),
            BinToPDFmap = cms.vstring("gaussPlusLinear")
        ),
        Glb_pt_eta_mcTrue = cms.PSet(
            EfficiencyCategoryAndState = cms.vstring("passingGlb","pass"),
            UnbinnedVariables = cms.vstring("mass"),
            BinnedVariables = cms.PSet(
                mcTrue = cms.vstring("true"),
                pt = PTBINS,
                eta = ETABINS
            ),
            BinToPDFmap = cms.vstring()
        ),
        Trk_pt_eta = cms.PSet(
            EfficiencyCategoryAndState = cms.vstring("passingTrkEx","pass"),
            UnbinnedVariables = cms.vstring("mass"),
            BinnedVariables = cms.PSet(
                pt = PTBINS,
                eta = ETABINS
            ),
            BinToPDFmap = cms.vstring("gaussPlusLinear")
        ),
        Trk_pt_eta_mcTrue = cms.PSet(
            EfficiencyCategoryAndState = cms.vstring("passingTrkEx","pass"),
            UnbinnedVariables = cms.vstring("mass"),
            BinnedVariables = cms.PSet(
                mcTrue = cms.vstring("true"),
                pt = PTBINS,
                eta = ETABINS
            ),
            BinToPDFmap = cms.vstring()
        ),
    )
)

process.MC_MuFromTk = Template.clone(
    InputFileNames = cms.vstring("tnpJPsi_MuonID_JPsi-0.084pb.root", "tnpJPsi_MuonID_ppMuX.root"),
    InputDirectoryName = cms.string("histoMuFromTk"),
    InputTreeName = cms.string("mcUnbias_tree"),
    OutputFileName = cms.string(FILEPREFIX+"MC_MuFromTk.root"),
    Efficiencies = cms.PSet(
        Glb_pt_eta = cms.PSet(
            EfficiencyCategoryAndState = cms.vstring("Glb","pass"),
            UnbinnedVariables = cms.vstring(),
            BinnedVariables = cms.PSet(
                pt = PTBINS,
                eta = ETABINS
            ),
            BinToPDFmap = cms.vstring()
        ),
        Trk_pt_eta = cms.PSet(
            EfficiencyCategoryAndState = cms.vstring("Trk","pass"),
            UnbinnedVariables = cms.vstring(),
            BinnedVariables = cms.PSet(
                pt = PTBINS,
                eta = ETABINS
            ),
            BinToPDFmap = cms.vstring()
        )
    )
)

process.MC_MuFromCal = process.MC_MuFromTk.clone(
    InputDirectoryName = cms.string("histoMuFromCal"),
    OutputFileName = cms.string(FILEPREFIX+"MC_MuFromCal.root"),
    Efficiencies = cms.PSet(
        Glb_pt_eta = cms.PSet(
            EfficiencyCategoryAndState = cms.vstring("passingGlb","pass"),
            UnbinnedVariables = cms.vstring(),
            BinnedVariables = cms.PSet(
                pt = PTBINS,
                eta = ETABINS
            ),
            BinToPDFmap = cms.vstring()
        ),
        Trk_pt_eta = cms.PSet(
            EfficiencyCategoryAndState = cms.vstring("passingTrkEx","pass"),
            UnbinnedVariables = cms.vstring(),
            BinnedVariables = cms.PSet(
                pt = PTBINS,
                eta = ETABINS
            ),
            BinToPDFmap = cms.vstring()
        )
    )
)

# -------------TRIGGER-----------

process.TnP_TriggerFromGlb = Template.clone(
    InputFileNames = cms.vstring("tnpJPsi_Trigger_JPsi-0.084pb.root", "tnpJPsi_Trigger_ppMuX.root"),
    InputDirectoryName = cms.string("histo"),
    InputTreeName = cms.string("fitter_tree"),
    OutputFileName = cms.string(FILEPREFIX+"TnP_TriggerFromGlb.root"),
    Efficiencies = cms.PSet()
)

process.TnP_TriggerFromTrk = process.TnP_TriggerFromGlb.clone(
    OutputFileName = cms.string(FILEPREFIX+"TnP_TriggerFromTrk.root"),
)

process.MC_TriggerFromGlb = process.TnP_TriggerFromGlb.clone(
    InputTreeName = cms.string("mcUnbias_tree"),
    OutputFileName = cms.string(FILEPREFIX+"MC_TriggerFromGlb.root"),
)

process.MC_TriggerFromTrk = process.MC_TriggerFromGlb.clone(
    OutputFileName = cms.string(FILEPREFIX+"MC_TriggerFromTrk.root"),
)

for trig in [ 'HLTMu3', 'L1DiMuOpen' ]:
    setattr( process.TnP_TriggerFromGlb.Efficiencies, trig+'pt_eta',
	cms.PSet(
            EfficiencyCategoryAndState = cms.vstring(trig,"pass"),
            UnbinnedVariables = cms.vstring("mass"),
            BinnedVariables = cms.PSet(
                Glb = cms.vstring("pass"),
                pt = PTBINS,
                eta = ETABINS
            ),
            BinToPDFmap = cms.vstring("gaussPlusLinear")
        )
    )
    setattr( process.TnP_TriggerFromGlb.Efficiencies, trig+'pt_eta_mcTrue',
        cms.PSet(
            EfficiencyCategoryAndState = cms.vstring(trig,"pass"),
            UnbinnedVariables = cms.vstring("mass"),
            BinnedVariables = cms.PSet(
                Glb = cms.vstring("pass"),
                mcTrue = cms.vstring("true"),
                pt = PTBINS,
                eta = ETABINS
            ),
            BinToPDFmap = cms.vstring()
        )
    )
    setattr( process.MC_TriggerFromGlb.Efficiencies, trig+'pt_eta',
        cms.PSet(
            EfficiencyCategoryAndState = cms.vstring(trig,"pass"),
            UnbinnedVariables = cms.vstring(),
            BinnedVariables = cms.PSet(
                Glb = cms.vstring("pass"),
                pt = PTBINS,
                eta = ETABINS
            ),
            BinToPDFmap = cms.vstring()
        )
    )
    setattr( process.TnP_TriggerFromTrk.Efficiencies, trig+'pt_eta',
	cms.PSet(
            EfficiencyCategoryAndState = cms.vstring(trig,"pass"),
            UnbinnedVariables = cms.vstring("mass"),
            BinnedVariables = cms.PSet(
                Trk = cms.vstring("pass"),
                pt = PTBINS,
                eta = ETABINS
            ),
            BinToPDFmap = cms.vstring("gaussPlusLinear")
        )
    )
    setattr( process.TnP_TriggerFromTrk.Efficiencies, trig+'pt_eta_mcTrue',
        cms.PSet(
            EfficiencyCategoryAndState = cms.vstring(trig,"pass"),
            UnbinnedVariables = cms.vstring("mass"),
            BinnedVariables = cms.PSet(
                Trk = cms.vstring("pass"),
                mcTrue = cms.vstring("true"),
                pt = PTBINS,
                eta = ETABINS
            ),
            BinToPDFmap = cms.vstring()
        )
    )
    setattr( process.MC_TriggerFromTrk.Efficiencies, trig+'pt_eta',
        cms.PSet(
            EfficiencyCategoryAndState = cms.vstring(trig,"pass"),
            UnbinnedVariables = cms.vstring(),
            BinnedVariables = cms.PSet(
                Trk = cms.vstring("pass"),
                pt = PTBINS,
                eta = ETABINS
            ),
            BinToPDFmap = cms.vstring()
        )
    )

process.p = cms.Path(
    process.TnP_MuFromTk
    + process.MC_MuFromTk
    + process.TnP_MuFromCal
    + process.MC_MuFromCal
    + process.TnP_TriggerFromGlb
    + process.MC_TriggerFromGlb
    + process.TnP_TriggerFromTrk
    + process.MC_TriggerFromTrk
)


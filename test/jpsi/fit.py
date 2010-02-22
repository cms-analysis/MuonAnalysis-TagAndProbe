import FWCore.ParameterSet.Config as cms
process = cms.Process("Fit")

## Note: mass ranges for the fit don't have to be the same as the mass ranges for the 
##       making of T&P pairs  in tp_from_skim.py
massRange = ( 2.8, 3.5 ); massRangeSta = (2, 4.5)

# Add your own files here
process.source = cms.Source("EmptySource")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )    

def defineEfficiencies(Module, ptBins=None, etaBins=None, passing="passing", model=None, mc=False):
    Module.Categories = cms.untracked.PSet()
    setattr(Module.Categories, passing, cms.untracked.vstring(passing, "dummy[pass=1,fail=0]"))
    if mc:
        Module.Categories.mcTrue = cms.untracked.vstring("MC true", "dummy[true=1,false=0]")
    # defines a set of efficiency calculations, what PDF to use for fitting and how to bin the data;
    # there will be a separate output directory for each calculation that includes a simultaneous fit, side band subtraction and counting. 
    Module.Efficiencies = cms.untracked.PSet(
        #the name of the parameter set becomes the name of the directory
        pt = cms.untracked.PSet(
            #specifies the efficiency of which category and state to measure 
            EfficiencyCategoryAndState = cms.untracked.vstring(passing,"pass"),
            #specifies what unbinned variables to include in the dataset, the mass is needed for the fit
            UnbinnedVariables = cms.untracked.vstring("mass"),
            #specifies the binning of parameters
            BinnedVariables = cms.untracked.PSet(pt = ptBins),
            #first string is the default followed by binRegExp - PDFname pairs
            BinToPDFmap = cms.untracked.vstring(model)
        ),
        eta = cms.untracked.PSet(
            #specifies the efficiency of which category and state to measure 
            EfficiencyCategoryAndState = cms.untracked.vstring(passing,"pass"),
            #specifies what unbinned variables to include in the dataset, the mass is needed for the fit
            UnbinnedVariables = cms.untracked.vstring("mass"),
            #specifies the binning of parameters
            BinnedVariables = cms.untracked.PSet(eta = etaBins),
            #first string is the default followed by binRegExp - PDFname pairs
            BinToPDFmap = cms.untracked.vstring(model)
        ),

        pt_eta = cms.untracked.PSet(
            EfficiencyCategoryAndState = cms.untracked.vstring(passing,"pass"),
            UnbinnedVariables = cms.untracked.vstring("mass"),
            BinnedVariables = cms.untracked.PSet(
                pt = ptBins,
                eta = etaBins,
            ),
            BinToPDFmap = cms.untracked.vstring(model)
        ),
    )
    if mc:
        Module.Efficiencies.pt_mcTrue = cms.untracked.PSet(
            EfficiencyCategoryAndState = cms.untracked.vstring(passing,"pass"),
            UnbinnedVariables = cms.untracked.vstring("mass"),
            BinnedVariables = cms.untracked.PSet(
                mcTrue = cms.untracked.vstring("true"),
                pt = ptBins,
            )
            #unspecified binToPDFmap means no fitting
        )
        Module.Efficiencies.eta_mcTrue = cms.untracked.PSet(
            EfficiencyCategoryAndState = cms.untracked.vstring(passing,"pass"),
            UnbinnedVariables = cms.untracked.vstring("mass"),
            BinnedVariables = cms.untracked.PSet(
                mcTrue = cms.untracked.vstring("true"),
                eta = etaBins, 
            )
            #unspecified binToPDFmap means no fitting
        )
        Module.Efficiencies.pt_eta_mcTrue = cms.untracked.PSet(
            EfficiencyCategoryAndState = cms.untracked.vstring(passing,"pass"),
            UnbinnedVariables = cms.untracked.vstring("mass"),
            BinnedVariables = cms.untracked.PSet(
                mcTrue = cms.untracked.vstring("true"),
                pt = ptBins,
                eta = etaBins,
            )
        )
    

RunFit  = cms.EDAnalyzer("TagProbeFitTreeAnalyzer",
    # IO parameters:
    InputFileName = cms.untracked.string("tnpJPsi.root"),
    InputTreeName = cms.untracked.string("fitter_tree"),
    #numbrer of CPUs to use for fitting
    NumCPU = cms.untracked.uint32(4),
    # specifies wether to save the RooWorkspace containing the data for each bin and
    # the pdf object with the initial and final state snapshots
    SaveWorkspace = cms.untracked.bool(False),

    # defines all the real variables of the probes available in the input tree and intended for use in the efficiencies
    Variables = cms.untracked.PSet(
        mass = cms.untracked.vstring("Tag-Probe Mass", str(massRange[0]), str(massRange[1]), "GeV/c^{2}"),
        pt = cms.untracked.vstring("Probe p_{T}", "0", "1000", "GeV/c"),
        eta = cms.untracked.vstring("Probe #eta", "-2.4", "2.4", "")
    ),

    # defines all the discrete variables of the probes available in the input tree and intended for use in the efficiency calculations
    Categories = cms.untracked.PSet(
        mcTrue = cms.untracked.vstring("MC true", "dummy[true=1,false=0]"),
        passing = cms.untracked.vstring("isMuon", "dummy[pass=1,fail=0]")
    ),

    # defines all the PDFs that will be available for the efficiency calculations; uses RooFit's "factory" syntax;
    # each pdf needs to define "signal", "backgroundPass", "backgroundFail" pdfs, "efficiency[0.9,0,1]" and "signalFractionInPassing[0.9]" are used for initial values  
    PDFs = cms.untracked.PSet(
        gaussPlusFloatExpo = cms.untracked.vstring(
            "Gaussian::signal(mass, mean[3.09,2.9,3.2], sigma[0.03,0.01,0.05])",
            "Exponential::backgroundPass(mass, slopePass[0,-10,10])",
            "Exponential::backgroundFail(mass, slopeFail[0,-10,10])",
            "efficiency[0.9,0,1]",
            "signalFractionInPassing[0.9]"
        ),
        gaussPlusChebychev3 = cms.untracked.vstring(
            "Gaussian::signal(mass, mean[3.09,2.9,3.2], sigma[0.15,0.05,1.0])",
            "Chebychev::backgroundPass(mass, {c1[0,-1,1], c2[0,-1,1], c3[0,-1,1]})",
            "Chebychev::backgroundFail(mass, {c1, c2, c3})",
            "efficiency[0.9,0,1]",
            "signalFractionInPassing[0.9]"
        )
    ),

    # defines a set of efficiency calculations, what PDF to use for fitting and how to bin the data;
    # there will be a separate output directory for each calculation that includes a simultaneous fit, side band subtraction and counting. 
 )


process.fitGlbFromTk = RunFit.clone(
    InputDirectoryName = cms.untracked.string('histoMuFromTk'),
    OutputFileName     = cms.untracked.string('fitGlbFromTk.root'),
)
defineEfficiencies(process.fitGlbFromTk,
    ptBins  = cms.untracked.vdouble( 1.5, 3, 4.5, 6, 10, 20),
    etaBins = cms.untracked.vdouble( -2.4, -1.3, -0.8, 0.8, 1.3, 2.4),
    passing = 'passingGlb',
    model = 'gaussPlusFloatExpo',
    mc = True
)

process.fitGlbFromCal = process.fitGlbFromTk.clone(
    InputDirectoryName = 'histoMuFromCal',
    OutputFileName     = 'fitGlbFromCal.root',
)
process.fitTrkExFromTk  = RunFit.clone(
    InputDirectoryName = cms.untracked.string('histoMuFromTk'),
    OutputFileName     = cms.untracked.string('fitTrkExFromTk.root'),
)
defineEfficiencies(process.fitTrkExFromTk,
    ptBins  = cms.untracked.vdouble( 1.5, 3, 4.5, 6, 10, 20),
    etaBins = cms.untracked.vdouble( -2.4, -1.3, -0.8, 0.8, 1.3, 2.4),
    passing = 'passingTrkEx',
    model = 'gaussPlusFloatExpo',
    mc = True
)
process.fitTrkExFromCal = process.fitTrkExFromTk.clone(
    InputDirectoryName = 'histoMuFromCal',
    OutputFileName     = 'fitTrkExFromCal.root' 
)

process.fitTkFromSta = RunFit.clone(
    InputDirectoryName = cms.untracked.string('histoTkFromSta'),
    OutputFileName     = cms.untracked.string('fitTkFromSta.root'),
)
defineEfficiencies(process.fitTkFromSta,
    ptBins  = cms.untracked.vdouble( 1.5, 3, 4.5, 6, 10, 20),
    etaBins = cms.untracked.vdouble( -2.4, -1.3, -0.8, 0.8, 1.3, 2.4),
    passing = 'passing',
    model = 'gaussPlusChebychev3',
    mc = True
)
process.fitTkFromSta.Variables.mass = cms.untracked.vstring("Tag-Probe Mass", str(massRangeSta[0]), str(massRangeSta[1]), "GeV/c^{2}")

for L in [ 'Mu3', 'Mu5', 'DoubleMu0', 'DoubleMu3', 'L1MuOpen' ]:
    label = "fit%sFromGlb" % (L,);
    setattr(process,
        label,
        RunFit.clone(
            InputDirectoryName = cms.untracked.string('histoHltFromGlb'),
            OutputFileName     = cms.untracked.string(label+'.root'),
        )
    )
    getattr(process,label).Variables.eta = cms.untracked.vstring("Probe #eta", "-2.1", "2.1", "")
    defineEfficiencies(getattr(process,label),
        ptBins  = cms.untracked.vdouble( 1.5, 3, 4.5, 6, 10, 20),
        etaBins = cms.untracked.vdouble( -2.1,-1.2,-0.7,0.0,0.7,1.2,2.1),
        passing = L,
        model = 'gaussPlusFloatExpo',
        mc = True
    )

process.fitness = cms.Path(
    process.fitGlbFromTk +
    process.fitGlbFromCal +
    process.fitTrkExFromTk +
    process.fitTrkExFromCal +
    process.fitTkFromSta +
    process.fitMu3FromGlb +
    process.fitMu5FromGlb +
    process.fitDoubleMu0FromGlb +
    process.fitDoubleMu3FromGlb +
    process.fitL1MuOpenFromGlb 
)

#process.TFileService = cms.Service("TFileService", fileName = cms.string("fitJPsi.root"))

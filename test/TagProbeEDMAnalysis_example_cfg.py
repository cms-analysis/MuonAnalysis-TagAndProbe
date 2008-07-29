import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")
# keep the logging output to a nice level
process.load("FWCore.MessageService.MessageLogger_cfi")

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:output/generic_tp_edm_ntuple_cmssw207_1.root')
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)
process.demo = cms.EDFilter("TagProbeEDMAnalysis",
    # Efficiency/Fitting variables
    CalculateEffSideBand = cms.untracked.bool(True),
    NameVar2 = cms.untracked.string('eta'),
    # Variables for sideband subtraction
    SBSPeak = cms.untracked.double(91.1876),
    # Variable Specifications for SB subtractions and Roofit
    # Choose binned or unbinned fitter ...
    # Note that the unbinned fit will not fit weighted data,
    # if you wish to use weights, use the binned fit.
    UnbinnedFit = cms.untracked.bool(False),
    # Variables and binning for the eff hists
    # Valid variables names for the eff binning are:
    # "pt","p","px","py","pz","e","et","eta" and "phi"
    # If omitted the defaults are var1 = pt and var2 = eta
    NameVar1 = cms.untracked.string('pt'),
    # This way of declaring the bin will overide any other
    Var1BinBoundaries = cms.untracked.vdouble(0.0, 30.0, 40.0, 50.0, 100.0),
    NumBkgPass = cms.untracked.vdouble(1000.0, -10.0, 10000.0),
    BkgBeta = cms.untracked.vdouble(0.001, 0.0, 0.1),
    NumBinsVar2 = cms.untracked.int32(1),
    # Background variables
    BkgAlpha = cms.untracked.vdouble(62.0, 50.0, 70.0),
    # Root file to eff histograms to
    FitFileName = cms.untracked.string('muon_eff_test.root'),
    BkgPeak = cms.untracked.vdouble(91.1876),
    # Binning for the plots ...
    XBins = cms.untracked.vuint32(100, 100, 100, 100),
    SignalWidth = cms.untracked.vdouble(2.8, 1.0, 4.0),
    MassLow = cms.untracked.double(50.0),
    # Mass window for fitting
    NumBinsMass = cms.untracked.int32(50),
    logY = cms.untracked.vuint32(1, 1, 1, 0),
    # Efficiency variables
    Efficiency = cms.untracked.vdouble(0.98, 0.0, 1.0),
    BkgGamma = cms.untracked.vdouble(0.05, 0.0, 0.1),
    # Request a 2D fit as well if desired ... will add 2D eff hists
    # binned in Pt and Var2
    Do2DFit = cms.untracked.bool(True),
    XMax = cms.untracked.vdouble(120.0, 120.0, 120.0, 3.14),
    outputFileNames = cms.untracked.vstring('Zmass.gif', 
        'Zmass_pass.gif', 
        'MCZmass.gif', 
        'Phi.gif'),
    CalculateEffFitter = cms.untracked.bool(True), ## Calculate and store effs from Roofit

    SignalSigma = cms.untracked.vdouble(1.1, 0.5, 4.0),
    # Make some plots of tree variables ...
    quantities = cms.untracked.vstring('TPmass', 
        'TPmass', 
        'MC23mass', 
        'TPProbephi'),
    Var2High = cms.untracked.double(2.5),
    MassHigh = cms.untracked.double(120.0),
    BifurGaussFrac = cms.untracked.vdouble(0.87),
    #untracked vdouble Var2BinBoundaries   = {-2.5,2.5}
    # Fitter variables - for the Roofit fitter
    # If you want the variable to float in the fit fill
    # three array elements {default, range_low, range_high}
    # If the variable should be fixed, fill one element {value}
    # Signal variables
    SignalMean = cms.untracked.vdouble(91.1876, 90.0, 92.0),
    XMin = cms.untracked.vdouble(50.0, 50.0, 50.0, -3.14),
    NumBkgFail = cms.untracked.vdouble(10.0, -10.0, 10000.0),
    conditions = cms.untracked.vstring('1', 
        'TPppass==1', 
        '1', 
        '1'),
    SignalWidthL = cms.untracked.vdouble(3.0, 0.0, 20.0),
    SignalWidthR = cms.untracked.vdouble(0.52, 0.0, 10.0),
    Var2Low = cms.untracked.double(-2.5),
    SBSStanDev = cms.untracked.double(2.0), ## SD from peak for subtraction

    CalculateEffTruth = cms.untracked.bool(True), ## Calculate and store true effs

    NumSignal = cms.untracked.vdouble(4000.0, -10.0, 100000.0)
)

process.p = cms.Path(process.demo)
process.MessageLogger.destinations = ['cout']


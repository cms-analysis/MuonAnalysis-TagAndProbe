import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")
# keep the logging output to a nice level
process.load("FWCore.MessageService.MessageLogger_cfi")

process.source = cms.Source("EmptySource")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1)
)
process.demo = cms.EDFilter("TagProbeAnalysis",
    # Efficiency/Fitting variables
    CalculateEffSideBand = cms.untracked.bool(True),
    # Variables for sideband subtraction
    SBSPeak = cms.untracked.double(91.1876),
    # Variable Specifications for SB subtractions and Roofit
    # Choose binned or unbinned fitter ...
    # Note that the unbinned fit will not fit weighted data,
    # if you wish to use weights, use the binned fit.
    UnbinnedFit = cms.untracked.bool(False),
    NumBkgPass = cms.untracked.vdouble(1000.0, -10.0, 10000.0),
    BkgBeta = cms.untracked.vdouble(0.001, 0.0, 0.1),
    # Background variables
    BkgAlpha = cms.untracked.vdouble(62.0, 50.0, 70.0),
    # Root file to eff histograms to
    FitFileName = cms.untracked.string('muon_eff_test.root'),
    BkgPeak = cms.untracked.vdouble(91.1876),
    # Binning for the plots ...
    XBins = cms.untracked.vuint32(100, 100, 100),
    SignalWidth = cms.untracked.vdouble(2.8, 1.0, 4.0),
    EtaLow = cms.untracked.double(-2.5),
    MassLow = cms.untracked.double(50.0),
    NumBinsEta = cms.untracked.int32(1),
    # Mass window for fitting
    NumBinsMass = cms.untracked.int32(50),
    logY = cms.untracked.vuint32(1, 1, 1),
    # Efficiency variables
    Efficiency = cms.untracked.vdouble(0.98, 0.0, 1.0),
    BkgGamma = cms.untracked.vdouble(0.05, 0.0, 0.1),
    XMax = cms.untracked.vdouble(120.0, 120.0, 120.0),
    outputFileNames = cms.untracked.vstring('Zmass.gif', 
        'Zmass_pass.gif', 
        'MCZmass.gif'),
    CalculateEffFitter = cms.untracked.bool(True), ## Calculate and store effs from Roofit

    PtLow = cms.untracked.double(0.0),
    EtaHigh = cms.untracked.double(2.5),
    PtHigh = cms.untracked.double(100.0),
    SignalSigma = cms.untracked.vdouble(1.1, 0.5, 4.0),
    # Make some plots of tree variables ...
    quantities = cms.untracked.vstring('TPmass', 
        'TPmass', 
        'MC23mass'),
    MassHigh = cms.untracked.double(120.0),
    BifurGaussFrac = cms.untracked.vdouble(0.87),
    # Fitter variables - for the Roofit fitter
    # If you want the variable to float in the fit fill
    # three array elements {default, range_low, range_high}
    # If the variable should be fixed, fill one element {value}
    # Signal variables
    SignalMean = cms.untracked.vdouble(91.1876, 90.0, 92.0),
    XMin = cms.untracked.vdouble(50.0, 50.0, 50.0),
    NumBkgFail = cms.untracked.vdouble(10.0, -10.0, 10000.0),
    conditions = cms.untracked.vstring('1', 
        'TPppass==1', 
        '1'),
    SignalWidthR = cms.untracked.vdouble(0.52, 0.0, 10.0),
    SignalWidthL = cms.untracked.vdouble(3.0, 0.0, 20.0),
    SBSStanDev = cms.untracked.double(2.0), ## SD from peak for subtraction

    # Enter the input files here .. you can have any number of files with any normalization.
    # If you have already combined signal and background with the right normalizations (or have data!!)
    # just enter one filename. If you want to normalize signal and background samples differently,
    # enter the file-set for each sample with the appropriate cross-section and desired
    # integrated luminosity. In this example we enter signal files and background files, 
    # and normalize the signal to the 2.18pb^-1 of background.
    #untracked vstring inputFileNames = {
    #  "output/ZmumuSignalOut/zmumu_reco_muon_eff_*.root",
    #  "output/FilteredQCDOut/mu_pair_tagprobe_reco_30_120_*.root"
    #}
    # Enter a weight for each file name in the form of lumi*cross-section
    # if the files are desired to be unweighted enter -1.0 as the value.
    # If no weights are entered, the default is unweighted.
    #untracked vdouble Luminosity      = {2.18,-1.0}
    #untracked vdouble CrossSection    = {829.5,-1.0}
    inputFileNames = cms.untracked.vstring('output/generic_tp_edm_ntuple_cmssw207_1.root'),
    # Eta and Pt binning for the eff hists
    NumBinsPt = cms.untracked.int32(1),
    CalculateEffTruth = cms.untracked.bool(True), ## Calculate and store true effs

    NumSignal = cms.untracked.vdouble(4000.0, -10.0, 100000.0)
)

process.p = cms.Path(process.demo)
process.MessageLogger.destinations = ['cout', 'cerr']


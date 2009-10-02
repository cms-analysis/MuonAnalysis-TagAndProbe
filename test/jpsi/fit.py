import FWCore.ParameterSet.Config as cms
process = cms.Process("Fit")

massRange = ( 2.8, 3.5 ); massRangeSta = (2, 5)

# Add your own files here
process.source = cms.Source("EmptySource")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )    

RunFit = cms.EDAnalyzer("TagProbeEDMAnalysis",  
      ## Efficiency/Fitting variables
      CalculateEffSideBand = cms.untracked.bool( True ), ## Calculate and store effs using SB
      CalculateEffFitter   = cms.untracked.bool( True ), ## Calculate and store effs from Roofit
      CalculateEffTruth    = cms.untracked.bool( True ), ## Calculate and store true effs

      ## Set mode to read from files ...
      Mode = cms.untracked.string("Read"),

      ReadFromFiles = cms.untracked.vstring("TO_BE_FILLED"),
      FitFileName   = cms.untracked.string("TO_BE_FILLED"),

      UnbinnedFit  = cms.untracked.bool( True ),
      Do2DFit      = cms.untracked.bool( True ),

      ## Mass window for fitting
      NumBinsMass         = cms.untracked.int32( 20 ), # matters only for the plots, the fit is unbinned
      MassLow             = cms.untracked.double( massRange[0] ),
      MassHigh            = cms.untracked.double( massRange[1] ),
    
      ## Variable Specifications for SB subtractions and Roofit
      NameVar1             = cms.untracked.string( "pt" ),
      Var1BinBoundaries   = cms.untracked.vdouble( 3, 4.5, 6, 10, 20),
      NameVar2             = cms.untracked.string( "eta" ),
      Var2BinBoundaries   = cms.untracked.vdouble( -2.4,-1.2,-0.7,0.0,0.7,1.2,2.4),

      GaussLineShape = cms.untracked.PSet(
	GaussMean        = cms.untracked.vdouble( 3.09, 2.9,  3.1  ),
	GaussSigma       = cms.untracked.vdouble( 0.03, 0.01, 0.05 ),
      ),

      ## Background variables
      CMSBkgLineShape = cms.untracked.PSet(
	CMSBkgAlpha           = cms.untracked.vdouble( 0 ), # fix these two to zero
	CMSBkgBeta            = cms.untracked.vdouble( 0 ), # so it's just an exp
	CMSBkgPeak            = cms.untracked.vdouble( 3.09 ),
	CMSBkgGamma           = cms.untracked.vdouble( 0, -10, 10 )
      ),

      ## Efficiency variables
      Efficiency        = cms.untracked.vdouble( 0.90,0.0,1.0 ),    
      NumSignal         = cms.untracked.vdouble( 4000.0,-10.0,30000.0 ),    
      NumBkgPass        = cms.untracked.vdouble( 4000.0,-10.0,10000.0 ),    
      NumBkgFail        = cms.untracked.vdouble( 1000.0,-10.0,7000.0  ),    

      ## Variables for sideband subtraction
      SBSPeak            = cms.untracked.double( 3.09 ),   ## Mass peak
      SBSStanDev         = cms.untracked.double( 2 ),      ## SD from peak for subtraction

      # All the following is useless now that we just read and fit
      # but it's required...
      # --------------------------------------------
      MCTruthParentId = cms.untracked.int32(443),
)

process.fitGlbFromTk = RunFit.clone(
    ReadFromFiles = [ 'histo_output_GlbFromTk.root' ],
    FitFileName   =     'fit_result_GlbFromTk.root'  ,
)

process.fitGlbFromCal = RunFit.clone(
    ReadFromFiles = [ 'histo_output_GlbFromCal.root' ],
    FitFileName   =     'fit_result_GlbFromCal.root'  ,
)

process.fitTkFromSta = RunFit.clone(
    ReadFromFiles = [ 'histo_output_TkFromSta.root' ],
    FitFileName   =     'fit_result_TkFromSta.root'  ,
    ## Override mass window and line shape for Sta muons
    MassLow  = cms.untracked.double( massRangeSta[0] ),
    MassHigh = cms.untracked.double( massRangeSta[1] ),
    GaussLineShape = cms.untracked.PSet(
        GaussMean  = cms.untracked.vdouble( 3.09, 2.9,  3.1 ),
        GaussSigma = cms.untracked.vdouble( 0.15, 0.05, 1.0 ),
    ),
)

process.fitHltFromGlb = RunFit.clone(
    ReadFromFiles = [ 'histo_output_HltFromGlb.root' ],
    FitFileName   =     'fit_result_HltFromGlb.root'  ,
    ## Override eta bins for trigger
    Var2BinBoundaries   = cms.untracked.vdouble( -2.1,-1.2,-0.7,0.0,0.7,1.2,2.1),
)

process.fitness = cms.Path(
    process.fitGlbFromTk +
    process.fitGlbFromCal +
    process.fitTkFromSta +
    process.fitHltFromGlb 
)

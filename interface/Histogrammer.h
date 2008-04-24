// system include files
#include <memory>
#include <vector>
#include <string>
#include <algorithm>
#include <iterator>
#include <iostream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

// ROOT

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TTree.h>

class Histogrammer : public edm::EDAnalyzer
{
   public:
      explicit Histogrammer(const edm::ParameterSet&);
      ~Histogrammer();

   private:
      virtual void beginJob(const edm::EventSetup&) ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      int CreateHistogram(TH1F& Histo, int i);
      int SaveHistogram(TH1F& Histo, std::string outFileName, Int_t LogY);

      void SideBandSubtraction(const TH1F& Total, TH1F& Result, Double_t Peak, Double_t SD);
      void ZllEffFitter( TTree &fitTree, std::string &fileName, std::string &bvar,
			 int bnbins, double blow, double nhigh );
      void CalculateEfficiencies();
      void SBSExample();

      // Histogram drawing input variables
      std::vector<std::string> fileNames;
      std::vector<std::string> quantities;  
      std::vector<std::string> conditions;  
      std::vector<std::string> outputFileNames;

      std::vector<unsigned int> XBins;
      std::vector<double> XMax;
      std::vector<double> XMin;
      std::vector<unsigned int> logY;

      std::vector<double> lumi_;
      std::vector<double> xsection_;
      std::vector<double> weight_;

      // Fitting/Efficiency calculation inputs
      bool calcEffsSB_;         // Calculate effs using SB subtraction for these TTrees?
      bool calcEffsFitter_;     // Calculate effs using fitter for these TTrees?

      std::string fitFileName_; // Name of the root file to write eff info to

      int massNbins_;           // Number of bins in the fit
      double massLow_;          // Lower bound for fit range
      double massHigh_;         // Upper bound for fit range
      
      int ptNbins_;             // Number of pt eff bins
      double ptLow_;            // Lower bound for pt eff range
      double ptHigh_;           // Upper bound for pt eff range

      int etaNbins_;            // Number of eta eff bins
      double etaLow_;           // Lower bound for eta eff range
      double etaHigh_;          // Upper bound for eta eff range

      std::vector<double> signalMean_;       // Fit mean
      std::vector<double> signalWidth_;      // Fit width
      std::vector<double> signalSigma_;      // Fit sigma
      std::vector<double> signalWidthL_;     // Fit left width
      std::vector<double> signalWidthR_;     // Fit right width

      double bifurGaussFrac_;                // Fraction of signal shape from bifur Gauss

      double bkgAlpha_;                      // Fit background shape alpha
      double bkgBeta_;                       // Fit background shape beta
      double bkgPeak_;                       // Fit background shape peak
      std::vector<double> bkgGamma_;         // Fit background shape gamma

      std::vector<double> efficiency_;       // Signal efficiency from fit
      std::vector<double> numSignal_;        // Signal events from fit
      std::vector<double> numBkgPass_;       // Background events passing from fit
      std::vector<double> numBkgFail_;       // Background events failing from fit


      TChain* fChain;
      TH1F* Histograms;
      int* NumEvents;
  
      unsigned int vectorSize;
      bool doAnalyze_;
};

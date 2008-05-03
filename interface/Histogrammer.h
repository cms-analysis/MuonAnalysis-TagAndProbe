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
      void ZllEffFitter( TTree *fitTree, std::string &fileName, std::string &bvar,
			 int bnbins, double blow, double nhigh );
      void ZllEffSBS( TTree *fitTree, std::string &fileName, std::string &bvar,
			 int bnbins, double blow, double nhigh );
      void CalculateEfficiencies();
      void CalculateMCTruthEfficiencies();

      std::vector<std::string> fileNames_;

      // Histogram drawing input variables
      std::vector<std::string> quantities_;  
      std::vector<std::string> conditions_;  
      std::vector<std::string> outputFileNames_;

      std::vector<unsigned int> XBins_;
      std::vector<double> XMax_;
      std::vector<double> XMin_;
      std::vector<unsigned int> logY_;

      std::vector<double> lumi_;
      std::vector<double> xsection_;
      std::vector<double> weight_;

      // Fitting/Efficiency calculation inputs
      int tagProbeType_;        // If more than one tag-probe type is stored select

      bool calcEffsSB_;         // Calculate effs using SB subtraction for these TTrees?
      bool calcEffsFitter_;     // Calculate effs using fitter for these TTrees?
      bool calcEffsTruth_;      // Calculate effs using MC truth for these TTrees

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

      std::vector<double> bifurGaussFrac_;   // Fraction of signal shape from bifur Gauss

      std::vector<double> bkgAlpha_;         // Fit background shape alpha
      std::vector<double> bkgBeta_;          // Fit background shape beta
      std::vector<double> bkgPeak_;          // Fit background shape peak
      std::vector<double> bkgGamma_;         // Fit background shape gamma

      std::vector<double> efficiency_;       // Signal efficiency from fit
      std::vector<double> numSignal_;        // Signal events from fit
      std::vector<double> numBkgPass_;       // Background events passing from fit
      std::vector<double> numBkgFail_;       // Background events failing from fit

      double SBSPeak_;                       // Sideband sub peak
      double SBSStanDev_;                    // Sideband sub standard deviation

      TChain* fChain_;
      TH1F* Histograms_;
      int* NumEvents_;
  
      unsigned int numQuantities_;
      bool doAnalyze_;

      // Tree variables ... could use EvtTree class here
      Int_t           nrtp_;
      Int_t           tp_type_[100];
      Int_t           tp_true_[100];
      Int_t           tp_ppass_[100];
      Float_t         tp_mass_[100]; 
      Float_t         tp_dpt_[100][2];
      Float_t         tp_deta_[100][2];
      
      TBranch        *b_nrtp_;   
      TBranch        *b_tp_true_;   
      TBranch        *b_tp_type_;   
      TBranch        *b_tp_ppass_;   
      TBranch        *b_tp_mass_;   
      TBranch        *b_tp_dpt_;   
      TBranch        *b_tp_deta_;   

      Int_t           ncnd_;
      Int_t           cnd_type_[100];
      Int_t           cnd_tag_[100];
      Int_t           cnd_pprobe_[100];
      Int_t           cnd_aprobe_[100];
      Float_t         cnd_pt_[100]; 
      Float_t         cnd_eta_[100];
      
      TBranch        *b_ncnd_;   
      TBranch        *b_cnd_type_;   
      TBranch        *b_cnd_tag_;   
      TBranch        *b_cnd_pprobe_;   
      TBranch        *b_cnd_aprobe_;   
      TBranch        *b_cnd_pt_;   
      TBranch        *b_cnd_eta_;   
};

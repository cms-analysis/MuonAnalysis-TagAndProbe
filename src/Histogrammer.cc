// -*- C++ -*-
//
// Package:    Histogrammer
// Class:      Histogrammer
// 
/**\class Histogrammer Histogrammer.cc MuonAnalysis/Histogrammer/src/Histogrammer.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  "Adam Hunt"
//         Created:  Sun Apr 20 10:35:25 CDT 2008
// $Id: Histogrammer.cc,v 1.7 2008/04/24 19:45:09 ahunt Exp $
//
//

#include "MuonAnalysis/TagAndProbe/interface/Histogrammer.h"
#include "MuonAnalysis/TagAndProbe/interface/RooCMSShapePdf.h"


// ROOT headers

#include <TArrow.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TGaxis.h>
#include <TGraphAsymmErrors.h>
#include <TIterator.h>
#include <TLatex.h>
#include <TString.h>
#include <TStyle.h>

// RooFit headers

#include <RooAbsData.h>
#include <RooAddPdf.h>
#include <RooBifurGauss.h>
#include <RooBreitWigner.h>
#include <RooCategory.h>
#include <RooCatType.h>
#include <RooCBShape.h>
#include <RooChi2Var.h>
#include <RooDataSet.h>
#include <RooDataHist.h>
#include <RooFitResult.h>
#include <RooGenericPdf.h>
#include <RooGlobalFunc.h>
#include <RooLandau.h>
#include <RooMinuit.h>
#include <RooNLLVar.h>
#include <RooPlot.h>
#include <RooRealVar.h>
#include <RooSimultaneous.h>
#include <RooTreeData.h>
#include <RooVoigtian.h>

using namespace std;
using namespace RooFit;

Histogrammer::Histogrammer(const edm::ParameterSet& iConfig)
{

   fileNames_       = iConfig.getUntrackedParameter< vector<string> >("inputFileNames");

   // Histogrammer variables
   quantities_      = iConfig.getUntrackedParameter< vector<string> >("quantities"); 
   conditions_      = iConfig.getUntrackedParameter< vector<string> >("conditions"); 
   outputFileNames_ = iConfig.getUntrackedParameter< vector<string> >("outputFileNames");
   XBins_           = iConfig.getUntrackedParameter< vector<unsigned int> >("XBins");
   XMin_            = iConfig.getUntrackedParameter< vector<double> >("XMin");
   XMax_            = iConfig.getUntrackedParameter< vector<double> >("XMax");
   logY_            = iConfig.getUntrackedParameter< vector<unsigned int> >("logY");
   
   // Normalization variables
   lumi_           = iConfig.getUntrackedParameter< vector<double> >("Luminosity");
   xsection_       = iConfig.getUntrackedParameter< vector<double> >("CrossSection");

   // Efficiency input variables
   tagProbeType_   = iConfig.getUntrackedParameter< int >("TagProbeType",0);

   calcEffsSB_     = iConfig.getUntrackedParameter< bool >("CalculateEffSideBand",false);
   calcEffsFitter_ = iConfig.getUntrackedParameter< bool >("CalculateEffFitter",false);
   calcEffsTruth_  = iConfig.getUntrackedParameter< bool >("CalculateEffTruth",false);

   fitFileName_    = iConfig.getUntrackedParameter< string >("FitFileName","fitfile.root");

   massNbins_      = iConfig.getUntrackedParameter< int >("NumBinsMass",20);
   massLow_        = iConfig.getUntrackedParameter< double >("MassLow",0.0);
   massHigh_       = iConfig.getUntrackedParameter< double >("MassHigh",100.0);

   ptNbins_        = iConfig.getUntrackedParameter< int >("NumBinsPt",20);
   ptLow_          = iConfig.getUntrackedParameter< double >("PtLow",0.0);
   ptHigh_         = iConfig.getUntrackedParameter< double >("PtHigh",100.0);

   etaNbins_       = iConfig.getUntrackedParameter< int >("NumBinsEta",20);
   etaLow_         = iConfig.getUntrackedParameter< double >("EtaLow",-2.4);
   etaHigh_        = iConfig.getUntrackedParameter< double >("EtaHigh",2.4);

   // Fitter variables
   vector<double> dSigM;
   dSigM.push_back(91.1876);
   dSigM.push_back(85.0);
   dSigM.push_back(95.0);
   signalMean_     = iConfig.getUntrackedParameter< vector<double> >("SignalMean",dSigM);
   vector<double> dSigW;
   dSigW.push_back(2.3);
   dSigW.push_back(1.0);
   dSigW.push_back(4.0);
   signalWidth_     = iConfig.getUntrackedParameter< vector<double> >("SignalWidth",dSigW);
   vector<double> dSigS;
   dSigS.push_back(1.5);
   dSigS.push_back(0.0);
   dSigS.push_back(4.0);
   signalSigma_     = iConfig.getUntrackedParameter< vector<double> >("SignalSigma",dSigS);
   vector<double> dSigWL;
   dSigWL.push_back(3.0);
   dSigWL.push_back(1.0);
   dSigWL.push_back(10.0);
   signalWidthL_    = iConfig.getUntrackedParameter< vector<double> >("SignalWidthL",dSigWL);
   vector<double> dSigWR;
   dSigWR.push_back(0.52);
   dSigWR.push_back(0.0);
   dSigWR.push_back(2.0);
   signalWidthR_    = iConfig.getUntrackedParameter< vector<double> >("SignalWidthR",dSigWR);
   
   vector<double> dBGF;
   dBGF.push_back(0.87);
   dBGF.push_back(0.0);
   dBGF.push_back(1.0);
   bifurGaussFrac_  = iConfig.getUntrackedParameter< vector<double> >("BifurGaussFrac",dBGF);

   vector<double> dBAl;
   dBAl.push_back(63.0);
   bkgAlpha_        = iConfig.getUntrackedParameter< vector<double> >("BkgAlpha",dBAl);
   vector<double> dBBt;
   dBBt.push_back(0.001);
   bkgBeta_         = iConfig.getUntrackedParameter< vector<double> >("BkgBeta",dBBt);
   vector<double> dBPk;
   dBPk.push_back(91.1876);
   bkgPeak_         = iConfig.getUntrackedParameter< vector<double> >("BkgPeak",dBPk);
   vector<double> dBGam;
   dBGam.push_back(0.08);
   dBGam.push_back(0.0);
   dBGam.push_back(1.0);
   bkgGamma_        = iConfig.getUntrackedParameter< vector<double> >("BkgGamma",dBGam);

   vector<double> dEff;
   dEff.push_back(0.98);
   dEff.push_back(0.0);
   dEff.push_back(1.1);
   efficiency_      = iConfig.getUntrackedParameter< vector<double> >("Efficiency",dEff);
   vector<double> dNSig;
   dNSig.push_back(1000.0);
   dNSig.push_back(-10.0);
   dNSig.push_back(1000000.0);
   numSignal_       = iConfig.getUntrackedParameter< vector<double> >("NumSignal",dNSig);
   vector<double> dNBPs;
   dNBPs.push_back(1000.0);
   dNBPs.push_back(-10.0);
   dNBPs.push_back(1000000.0);
   numBkgPass_      = iConfig.getUntrackedParameter< vector<double> >("NumBkgPass",dNBPs);
   vector<double> dNBFl;
   dNBFl.push_back(1000.0);
   dNBFl.push_back(-10.0);
   dNBFl.push_back(1000000.0);
   numBkgFail_      = iConfig.getUntrackedParameter< vector<double> >("NumBkgFail",dNBFl);

   
   // Allocate space for histograms
   numQuantities_ = quantities_.size();
   Histograms_ = new TH1F[numQuantities_];

   // Chain files together
   std::string tempString;
   fChain_ = new TChain("evttree");
   NumEvents_ = new int[fileNames_.size()];

   vector<string>::iterator it=fileNames_.begin();
   for( int i = 0; it < fileNames_.end(); it++, i++)
   {
      tempString = *it;
      fChain_->Add(tempString.c_str());
      NumEvents_[i] = fChain_->GetEntries();
   }

   fChain_->SetBranchStatus("*",0);
   if( calcEffsFitter_ )
   {
      fChain_->SetBranchAddress("nrtp", &nrtp_, &b_nrtp_);
      fChain_->SetBranchAddress("tp_true", tp_true_, &b_tp_true_);
      fChain_->SetBranchAddress("tp_type", tp_type_, &b_tp_type_);
      fChain_->SetBranchAddress("tp_ppass", tp_ppass_, &b_tp_ppass_);
      fChain_->SetBranchAddress("tp_mass", tp_mass_, &b_tp_mass_);
      fChain_->SetBranchAddress("tp_dpt", tp_dpt_, &b_tp_dpt_);
      fChain_->SetBranchAddress("tp_deta", tp_deta_, &b_tp_deta_);
      
      fChain_->SetBranchStatus("nrtp",1);
      fChain_->SetBranchStatus("tp_type",1);
      fChain_->SetBranchStatus("tp_ppass",1);
      fChain_->SetBranchStatus("tp_mass",1);
      fChain_->SetBranchStatus("tp_dpt",1);
      fChain_->SetBranchStatus("tp_deta",1);
   }
   if( calcEffsTruth_ )
   {
      fChain_->SetBranchAddress("ncnd", &ncnd_, &b_ncnd_);
      fChain_->SetBranchAddress("cnd_tag", cnd_tag_, &b_cnd_tag_);
      fChain_->SetBranchAddress("cnd_type", cnd_type_, &b_cnd_type_);
      fChain_->SetBranchAddress("cnd_pprobe", cnd_pprobe_, &b_cnd_pprobe_);
      fChain_->SetBranchAddress("cnd_aprobe", cnd_aprobe_, &b_cnd_aprobe_);
      fChain_->SetBranchAddress("cnd_pt", cnd_pt_, &b_cnd_pt_);
      fChain_->SetBranchAddress("cnd_eta", cnd_eta_, &b_cnd_eta_);

      fChain_->SetBranchStatus("ncnd",1);
      fChain_->SetBranchStatus("cnd_type",1);
      fChain_->SetBranchStatus("cnd_pprobe",1);
      fChain_->SetBranchStatus("cnd_aprobe",1);
      fChain_->SetBranchStatus("cnd_pt",1);
      fChain_->SetBranchStatus("cnd_eta",1);
   }

   // Verify correct use of cfg
   doAnalyze_ = true;
   if(outputFileNames_.size() != numQuantities_){
      doAnalyze_ = false;
      cout << "outputFileNames is not the same size as quantities" << endl;    
   }else if(conditions_.size() != numQuantities_){
      doAnalyze_ = false;
      cout << "conditions is not the same size as quantities" << endl;    
   }else if(XBins_.size() != numQuantities_){
      doAnalyze_ = false;
      cout << "XBins is not the same size as quantities" << endl;    
   }else if(XMax_.size() != numQuantities_){
      doAnalyze_ = false;
      cout << "XMax is not the same size as quantities" << endl;    
   }else if(XMin_.size() != numQuantities_){
      doAnalyze_ = false;
      cout << "XMin is not the same size as quantities" << endl;    
   }else if(logY_.size() != numQuantities_){
      doAnalyze_ = false;
      cout << "logY is not the same size as quantities" << endl;    
   }

   // Set default lumi values to 10 pb^-1
   if( lumi_.size() != fileNames_.size() )
   {
      if( lumi_.size() > fileNames_.size() ) lumi_.clear();
      for(unsigned int i=lumi_.size(); i<fileNames_.size(); ++i ) lumi_.push_back(10.0);
   }
   // Set default cross-section values to 1000 pb
   if( xsection_.size() != fileNames_.size() )
   {
      if( xsection_.size() > fileNames_.size() ) xsection_.clear();
      for(unsigned int i=xsection_.size(); i<fileNames_.size(); ++i ) xsection_.push_back(1000.0);
   }

   // Set the raw weights ... i.e. the number of events to
   // scale to ...
   for(unsigned int i=0; i<lumi_.size(); ++i )
   {
      if( lumi_[i]>0.0 && xsection_[i]>0.0 ) weight_.push_back(lumi_[i]*xsection_[i]);
      else                                   weight_.push_back(1.0);
   }

}

Histogrammer::~Histogrammer()
{
  if (fChain_){
    delete fChain_->GetCurrentFile();
    fChain_ = 0;
  }
  if(Histograms_){
    delete [] Histograms_; 
    Histograms_ = 0;
  }
  if(NumEvents_){
    delete [] NumEvents_;
    NumEvents_ = 0;
  }
}

// ------------ method called to for each event  ------------
void
Histogrammer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  
   if( doAnalyze_ )
   {
      for(unsigned int i = 0; i < numQuantities_; i++)
      {
	 CreateHistogram(Histograms_[i], i);
	 SaveHistogram(Histograms_[i], outputFileNames_[i], logY_[i]);
      }
   }

   CalculateEfficiencies();
   //SBSExample();
}

// ******* Create the user requested histograms ******** //
int Histogrammer::CreateHistogram(TH1F& Histo, int i)
{
  
   std::stringstream tempsstream;
   std::stringstream HistoName;

   // initialize
   HistoName.str(std::string());
   HistoName << "Histo" << i;

   tempsstream.str(std::string());    
   tempsstream << quantities_[i] << " >>  " << HistoName.str();  

   if(fChain_->FindBranch(quantities_[i].substr(0,quantities_[i].find("[")).c_str())){ 
      Histo =  TH1F(HistoName.str().c_str(), "", XBins_[i], XMin_[i], XMax_[i]);
      fChain_->Draw(tempsstream.str().c_str(), conditions_[i].c_str());
   }else{
      std::cout << "Branch does not exist: " << quantities_[i] << std::endl;
   }

   return 0;
}

// ***************************************************** //

// ********* Save the user requested histograms ******** //
int Histogrammer::SaveHistogram(TH1F& Histo, std::string outFileName, Int_t LogY = 0)
{
  
   TCanvas* c1 = new TCanvas("c1","c1",700,500);
   c1->GetPad(0)->SetTicks(1,1);
   c1->SetLogy(LogY);
  
   Histo.Draw();
  
   c1->SaveAs(outFileName.c_str());
  
   delete c1;

   return 0;
}
// ***************************************************** //


// ********* Do sideband subtraction on the requested histogram ********* //
void Histogrammer::SideBandSubtraction( const TH1F& Total, TH1F& Result, 
					Double_t Peak, Double_t SD)
{
   // Total Means signal plus background

   const Double_t BinWidth  = Total.GetXaxis()->GetBinWidth(1);
   const Int_t nbins = Total.GetNbinsX();
   const Double_t xmin = Total.GetXaxis()->GetXmin();

   const Int_t PeakBin = (Int_t)((Peak - xmin)/BinWidth + 1); // Peak
   const Int_t SDBin = (Int_t)(SD/BinWidth); // Standard deviation
   const Int_t I = 3*SDBin; // Interval
   const Int_t D = 10*SDBin;  // Distance from peak

   const Double_t IntegralRight = Total.Integral(PeakBin + D, PeakBin + D + I);
   const Double_t IntegralLeft = Total.Integral(PeakBin - D - I, PeakBin - D);

   Double_t SubValue = 0.0;
   Double_t NewValue = 0.0;

   for(Int_t bin = 1; bin < (nbins + 1); bin++){
      SubValue = ((IntegralRight - IntegralLeft)/(2*D+I)*(bin - PeakBin - D - I/2.0) + IntegralRight)/I;
      if(SubValue < 0)
	 SubValue = 0;

      NewValue = Total.GetBinContent(bin)-SubValue;
      if(NewValue > 0){
	 Result.SetBinContent(bin, NewValue);
      }
   }
   Result.SetEntries(Result.Integral("width"));
}
// ********************************************************************** //

// ************ Sideband Example ************ //
void Histogrammer::SBSExample(){

  const Int_t SBSXBins = 200;
  const Int_t SBSXMin = 0;
  const Int_t SBSXMax = 200;

   TH1F TPMassHistoSig("TPMassHistoSig", "", SBSXBins, SBSXMin, SBSXMax);
   TH1F TPMassHistoBG("TPMassHistoBG",   "", SBSXBins, SBSXMin, SBSXMax);
   TH1F TPMassHistoTot("TPMassHistoTot", "", SBSXBins, SBSXMin, SBSXMax);
   TH1F TPMassHistoSBS("TPMassHistoSBS", "", SBSXBins, SBSXMin, SBSXMax);

   fChain_->Draw("tp_mass >> TPMassHistoSig","", "", NumEvents_[0], 0);
   fChain_->Draw("tp_mass >> TPMassHistoBG", "", "", NumEvents_[1] - NumEvents_[0], NumEvents_[0] + 1);
  
   const Double_t SigScale = weight_[0]/NumEvents_[0];

   std::cout << "SigScale: " << SigScale << std::endl;

   TPMassHistoTot.Add(&TPMassHistoSig, &TPMassHistoBG, SigScale, 1);
   SideBandSubtraction(TPMassHistoTot, TPMassHistoSBS, 90, 2);

   SaveHistogram(TPMassHistoSig, "TPMassHistoSig.png", 1);
   SaveHistogram(TPMassHistoBG,  "TPMassHistoBG.png",  1);
   SaveHistogram(TPMassHistoTot, "TPMassHistoTot.png", 1);
   SaveHistogram(TPMassHistoSBS, "TPMassHistoSBS.png", 1);

}
// ******************************************** //


// ********** Z -> l+l- Fitter ********** //
void Histogrammer::ZllEffFitter( TTree *fitTree, string &fileName, string &bvar,
				 int bnbins, double blow, double bhigh )
{
   TFile outRootFile(fileName.c_str(),"UPDATE");
   outRootFile.cd();

   //return;
   cout << "Here in Zll fitter" << endl;
   
   string hname = "heff_" + bvar;
   string htitle = "Efficiency vs " + bvar;
   TH1F effhist(hname.c_str(),htitle.c_str(),bnbins,blow,bhigh);

   double bwidth = (bhigh-blow)/(double)bnbins;

   for( int bin=0; bin<bnbins; ++bin )
   {

      // The fit variable - lepton invariant mass
      RooRealVar Mass("Mass","Invariant Di-Lepton Mass", massLow_, massHigh_, "GeV/c^{2}");
      Mass.setBins(massNbins_);

      // The binning variable
      string bunits = "";
      double lowEdge = blow + (double)bin*bwidth;
      double highEdge = lowEdge + bwidth;
      if( bvar == "Pt" ) bunits = "GeV";
      RooRealVar Pt(bvar.c_str(),bvar.c_str(),lowEdge,highEdge,bunits.c_str());

      // The weighting
      RooRealVar Weight("Weight","Weight",0.0,10000.0);

      // Make the category variable that defines the two fits,
      // namely whether the probe passes or fails the eff criteria.
      RooCategory ProbePass("ProbePass","sample");
      ProbePass.defineType("pass",1);
      ProbePass.defineType("fail",0);  

      cout << "Made fit variables" << endl;

      // Add the TTree as our data set ... with the weight in case 
      // we are using chained MC
      //RooDataSet* data = new RooDataSet("fitData","fitData",fitTree,
      //				RooArgSet(ProbePass,Mass,Pt,Weight),"","");
      // Above command doesn't work in root 5.18 (lovely) so we have this
      // silly workaround with TChain for now
      RooDataSet* data = new RooDataSet("fitData","fitData",fitTree,
					RooArgSet(ProbePass,Mass,Pt,Weight));


      //data->get()->Print();
      data->setWeightVar("Weight");
      data->get()->Print();

      cout << "Made dataset" << endl;

      RooDataHist *bdata = new RooDataHist("bdata","Binned Data",
					   RooArgList(Mass,ProbePass),*data);
 
      // ********** Construct signal shape PDF ********** //

      // Signal PDF variables
      RooRealVar signalMean("signalMean","signalMean",signalMean_[0]);
      RooRealVar signalWidth("signalWidth","signalWidth",signalWidth_[0]);
      RooRealVar signalSigma("signalSigma","signalSigma",signalSigma_[0]);
      RooRealVar signalWidthL("signalWidthL","signalWidthL",signalWidthL_[0]);
      RooRealVar signalWidthR("signalWidthR","signalWidthR",signalWidthR_[0]);

      // If the user has set a range, make the variable float
      if( signalMean_.size() == 3 )
      {
	 signalMean.setRange(signalMean_[1],signalMean_[2]);
	 signalMean.setConstant(false);
      }
      if( signalWidth_.size() == 3 )
      {
	 signalWidth.setRange(signalWidth_[1],signalWidth_[2]);
	 signalWidth.setConstant(false);
      }
      if( signalSigma_.size() == 3 )
      {
	 signalSigma.setRange(signalSigma_[1],signalSigma_[2]);
	 signalSigma.setConstant(false);
      }
      if( signalWidthL_.size() == 3 )
      {
	 signalWidthL.setRange(signalWidthL_[1],signalWidthL_[2]);
	 signalWidthL.setConstant(false);
      }
      if( signalWidthR_.size() == 3 )
      {
	 signalWidthR.setRange(signalWidthR_[1],signalWidthR_[2]);
	 signalWidthR.setConstant(false);
      }
  
      // Voigtian
      RooVoigtian signalVoigtPdf("signalVoigtPdf", "signalVoigtPdf", 
				 Mass, signalMean, signalWidth, signalSigma);

      // Bifurcated Gaussian
      RooBifurGauss signalGaussBifurPdf("signalGaussBifurPdf", "signalGaussBifurPdf", 
					Mass, signalMean, signalWidthL, signalWidthR);

      // Bifurcated Gaussian fraction
      RooRealVar bifurGaussFrac("bifurGaussFrac","bifurGaussFrac",bifurGaussFrac_[0]);
      if( bifurGaussFrac_.size() == 3 )
      {
	 bifurGaussFrac.setRange(bifurGaussFrac_[1],bifurGaussFrac_[2]);
	 bifurGaussFrac.setConstant(false);
      } 

      // The total signal PDF
      RooAddPdf  signalShapePdf("signalShapePdf", "signalShapePdf",
				signalVoigtPdf,signalGaussBifurPdf,bifurGaussFrac);

      // ********** Construct background shape PDF ********** //

      // Background PDF variables
      RooRealVar bkgAlpha("bkgAlpha","bkgAlpha",bkgAlpha_[0]);
      RooRealVar bkgBeta("bkgBeta","bkgBeta",bkgBeta_[0]);
      RooRealVar bkgGamma("bkgGamma","bkgGamma",bkgGamma_[0]);
      RooRealVar bkgPeak("bkgPeak","bkgPeak",bkgPeak_[0]);

      // If the user has specified a range, let the bkg shape 
      // variables float in the fit
      if( bkgAlpha_.size() == 3 )
      {
	 bkgAlpha.setRange(bkgAlpha_[1],bkgAlpha_[2]);
	 bkgAlpha.setConstant(false);
      }
      if( bkgBeta_.size() == 3 )
      {
	 bkgBeta.setRange(bkgBeta_[1],bkgBeta_[2]);
	 bkgBeta.setConstant(false);
      }
      if( bkgGamma_.size() == 3 )
      {
	 bkgGamma.setRange(bkgGamma_[1],bkgGamma_[2]);
	 bkgGamma.setConstant(false);
      }
      if( bkgPeak_.size() == 3 )
      {
	 bkgPeak.setRange(bkgPeak_[1],bkgPeak_[2]);
	 bkgPeak.setConstant(false);
      }

      // CMS Background shape
      RooCMSShapePdf bkgShapePdf("bkgShapePdf","bkgShapePdf", 
				 Mass,bkgAlpha,bkgBeta,bkgGamma,bkgPeak);

      // Now define some efficiency/yield variables  
      RooRealVar efficiency("efficiency","efficiency",efficiency_[0]);
      RooRealVar numSignal("numSignal","numSignal",numSignal_[0]);
      RooRealVar numBkgPass("numBkgPass","numBkgPass",numBkgPass_[0]);
      RooRealVar numBkgFail("numBkgFail","numBkgFail",numBkgFail_[0]);

      // If ranges are specifed these are floating variables
      if( efficiency_.size() == 3 )
      {
	 efficiency.setRange(efficiency_[1],efficiency_[2]);
	 efficiency.setConstant(false);
      }
      if( numSignal_.size() == 3 )
      {
	 numSignal.setRange(numSignal_[1],numSignal_[2]);
	 numSignal.setConstant(false);
      }
      if( numBkgPass_.size() == 3 )
      {
	 numBkgPass.setRange(numBkgPass_[1],numBkgPass_[2]);
	 numBkgPass.setConstant(false);
      }
      if( numBkgFail_.size() == 3 )
      {
	 numBkgFail.setRange(numBkgFail_[1],numBkgFail_[2]);
	 numBkgFail.setConstant(false);
      }
      

      RooFormulaVar numSigPass("numSigPass","numSignal*efficiency", 
			       RooArgList(numSignal,efficiency) );
      RooFormulaVar numSigFail("numSigFail","numSignal*(1.0 - efficiency)", 
			       RooArgList(numSignal,efficiency) );

      RooArgList componentspass(signalShapePdf,bkgShapePdf);
      RooArgList yieldspass(numSigPass, numBkgPass);
      RooArgList componentsfail(signalShapePdf,bkgShapePdf);
      RooArgList yieldsfail(numSigFail, numBkgFail);	  

      RooAddPdf sumpass("sumpass","fixed extended sum pdf",componentspass,yieldspass);
      RooAddPdf sumfail("sumfail","fixed extended sum pdf",componentsfail, yieldsfail);
  
      // The total simultaneous fit ...
      RooSimultaneous totalPdf("totalPdf","totalPdf",ProbePass);
      ProbePass.setLabel("pass");
      totalPdf.addPdf(sumpass,ProbePass.getLabel());
      totalPdf.Print();
      ProbePass.setLabel("fail");
      totalPdf.addPdf(sumfail,ProbePass.getLabel());
      totalPdf.Print();

      // Count the number of passing and failing probes in the region
      // making sure we have enough to fit ...
      cout << "About to count the number of events" << endl;
      int npassR = (int)data->sumEntries("ProbePass==1");
      int nfailR = (int)data->sumEntries("ProbePass==0");
      cout << "Num pass " << npassR << endl;
      cout << "Num fail " << nfailR << endl;

      RooAbsCategoryLValue& simCat = (RooAbsCategoryLValue&) totalPdf.indexCat();
   
      TList* dsetList = const_cast<RooAbsData*>((RooAbsData*)data)->split(simCat);
      RooCatType* type;
      TIterator* catIter = simCat.typeIterator();
      while( (type=(RooCatType*)catIter->Next()) )
      {
	 // Retrieve the PDF for this simCat state
	 RooAbsPdf* pdf =  totalPdf.getPdf(type->GetName());
	 RooAbsData* dset = (RooAbsData*) dsetList->FindObject(type->GetName());

	 if (pdf && dset && dset->numEntries(kTRUE)!=0.) 
	 {               
	    cout << "GOF Entries " << dset->numEntries() << " " 
		 << type->GetName() << std::endl;
	    if( (string)type->GetName() == "pass" ) 
	    {
	       npassR = dset->numEntries(); 
	       cout << "Pass " << npassR << endl; 
	    }
	    else if( (string)type->GetName() == "fail" ) 
	    {
	       nfailR = dset->numEntries();
	       cout << "Fail " << nfailR << endl; 
	    }
	 }
      }
      // End the pass fail counting.

      // Return if there's nothing to fit
      if( npassR==0 && nfailR==0 ) return;

      cout << "**** About to start the fitter ****" << endl;

      // ********* Do the Actual Fit ********** //  
      RooFitResult *fitResult = 0;
      //RooNLLVar nll("nll","nll",totalPdf,*bdata,kTRUE);
      //RooMinuit m(nll);
      RooChi2Var chi2("chi2","chi2",totalPdf,*bdata,
		      DataError(RooAbsData::SumW2),Extended(kTRUE));
      RooMinuit m(chi2);
      m.setErrorLevel(0.5); // <<< HERE
      m.setStrategy(2);
      m.hesse();
      m.migrad();
      m.hesse();
      m.minos();
      fitResult = m.save();

      fitResult->Print("v");

      std::cout << "Signal yield: " << numSignal.getVal() << " +- "
		<< numSignal.getError() << " + " << numSignal.getAsymErrorHi()
		<<" - "<< numSignal.getAsymErrorLo() << std::endl;
      std::cout << "Efficiency: "<< efficiency.getVal() << " +- "
		<< efficiency.getError() << " + " << efficiency.getAsymErrorHi()
		<<" + "<< efficiency.getAsymErrorLo() << std::endl;

      // Fill the efficiency hist
      effhist.SetBinContent(bin+1,efficiency.getVal());
      effhist.SetBinError(bin+1,efficiency.getError());

      // ********** Make and save Canvas for the plots ********** //
      outRootFile.cd();

      int font_num = 42;
      double font_size = 0.05;

      TStyle fitStyle("fitStyle","Style for Fit Plots");
      fitStyle.Reset("Plain");
      fitStyle.SetFillColor(10);
      fitStyle.SetTitleFillColor(10);
      fitStyle.SetTitleStyle(0000);
      fitStyle.SetStatColor(10);
      fitStyle.SetErrorX(0);
      fitStyle.SetEndErrorSize(10);
      fitStyle.SetPadBorderMode(0);
      fitStyle.SetFrameBorderMode(0);
      //fitStyle.SetOptTitle(0);

      fitStyle.SetTitleFont(font_num);
      fitStyle.SetTitleFontSize(font_size);
      fitStyle.SetTitleFont(font_num, "XYZ");
      fitStyle.SetTitleSize(font_size, "XYZ");
      fitStyle.SetTitleXOffset(0.9);
      fitStyle.SetTitleYOffset(1.05);
      fitStyle.SetLabelFont(font_num, "XYZ");
      fitStyle.SetLabelOffset(0.007, "XYZ");
      fitStyle.SetLabelSize(font_size, "XYZ");
      fitStyle.cd();

      ostringstream oss;
      oss << bin;
      string cname = "c_" + bvar + "_" + oss.str();
      TCanvas *c = new TCanvas(cname.c_str(),"Sum over Modes, Signal Region",1000,1500);
      c->Divide(1,2);
      c->cd(1);
      c->SetFillColor(10);
      TPad *lhs = (TPad*)gPad;
      lhs->Divide(2,1);
      lhs->cd(1);

      RooPlot* frame1 = Mass.frame();
      frame1->SetTitle("Passing Tag-Probes");
      frame1->SetName("pass");
      data->plotOn(frame1,Cut("ProbePass==1"));
      ProbePass.setLabel("pass");
      totalPdf.plotOn(frame1,Slice(ProbePass),ProjWData(Mass,*data));
      totalPdf.plotOn(frame1,Slice(ProbePass),Components(bkgShapePdf),
		      LineStyle(kDashed),ProjWData(Mass,*data));
      frame1->Draw("e0");

      outRootFile.cd();
      //frame1->Write();

      lhs->cd(2);
      RooPlot* frame2 = Mass.frame();
      frame2->SetTitle("Failing Tag-Probes");
      frame2->SetName("fail");
      data->plotOn(frame2,Cut("ProbePass==0"));
      ProbePass.setLabel("fail");
      totalPdf.plotOn(frame2,Slice(ProbePass),ProjWData(Mass,*data));
      totalPdf.plotOn(frame2,Slice(ProbePass),Components(bkgShapePdf),
		      LineStyle(kDashed),ProjWData(Mass,*data));
      frame2->Draw("e0");

      outRootFile.cd();
      //frame2->Write();

      c->cd(2);
      RooPlot* frame3 = Mass.frame();
      frame3->SetTitle("All Tag-Probes");
      frame3->SetName("total");
      data->plotOn(frame3);
      totalPdf.plotOn(frame3,ProjWData(Mass,*data));
      totalPdf.plotOn(frame3,Components(bkgShapePdf),
		      LineStyle(kDashed),ProjWData(Mass,*data));
      //frame3->Write();
      totalPdf.paramOn(frame3);
      frame3->Draw("e0");

      outRootFile.cd();
		
      std::cout << " Stop 0 " << std::endl;
      //c->Print("FitPlot.eps");
      outRootFile.cd();
      c->Write();

      std::cout << " Stop 1 " << std::endl;

      delete data;
      delete bdata;
   }

   outRootFile.cd();
   effhist.Write();

   //outRootFile.Write();
   outRootFile.Close();

   return;
}
// ************************************** //

// ********** Get the efficiency from this TTree ********** //
void Histogrammer::CalculateEfficiencies()
{
   if( calcEffsTruth_ ) CalculateMCTruthEfficiencies();

   if( calcEffsFitter_ )
   {
      // Loop over the number of different types of 
      // efficiency measurement in the input tree
      // Make a simple tree for fitting, and then
      // call the fitter.

      // Make the simple fit tree
      int    ProbePass;
      double Mass;
      double Pt;
      double Eta;
      double Weight;

      TTree *fitTree = new TTree("fitter_tree","Tree For Fitting",1);
      fitTree->Branch("ProbePass",&ProbePass,"ProbePass/I");
      fitTree->Branch("Mass",     &Mass,     "Mass/D");
      fitTree->Branch("Pt",       &Pt,       "Pt/D");
      fitTree->Branch("Eta",      &Eta,      "Eta/D");
      fitTree->Branch("Weight",   &Weight,   "Weight/D");

      int nFile = 0;
      Weight = 1.0;
      if( NumEvents_[nFile] > 0 ) Weight = weight_[nFile]/(double)NumEvents_[nFile];
      cout << "Filling fit tree with weight " << Weight << endl;
      for( int i=0; i<fChain_->GetEntries(); ++i )
      {
	 fChain_->GetEntry(i);

	 // If we pass a boundary, change the weight
	 if( nFile<(int)fileNames_.size() && i>=NumEvents_[nFile] )
	 {
	    ++nFile;
	    if( lumi_[nFile] > 0.0 )
	    {
	       if( NumEvents_[nFile] > 0 ) Weight = weight_[nFile]/(double)NumEvents_[nFile];
	    }
	    else
	    {
	       Weight = weight_[nFile];
	    }
	    cout << "Filling fit tree with weight " << Weight << endl;
	 }

	 for( int n=0; n<nrtp_; ++n )
	 {
	    if( tp_type_[n] != tagProbeType_ ) continue;

	    ProbePass = tp_ppass_[n];
	    Mass      = (double)tp_mass_[n];
	    Pt        = (double)tp_dpt_[n][1];
	    Eta       = (double)tp_deta_[n][1];

	    fitTree->Fill();
	 }
      }

      // Write the new TTree to the file for storage
      string fmode = "UPDATE";
      if( !calcEffsTruth_ ) fmode = "RECREATE";
      TFile outRootFile(fitFileName_.c_str(),fmode.c_str());
      outRootFile.cd();
      cout << "Writing tree" << endl; 
      fitTree->Write(); 
      outRootFile.Close();

      // We have filled the simple tree ... call the fitter
      string binnedVar = "Pt";
      ZllEffFitter( fitTree, fitFileName_, binnedVar, ptNbins_, ptLow_, ptHigh_ );
      binnedVar = "Eta";
      ZllEffFitter( fitTree, fitFileName_, binnedVar, etaNbins_, etaLow_, etaHigh_ );

   }

   if( calcEffsSB_ )
   {
      // Need to fill this function
   }

   return;
}
// ******************************************************** //

// ********** Get the true efficiency from this TTree ********** //
void Histogrammer::CalculateMCTruthEfficiencies()
{
   // Loop over the number of different types of 
   // efficiency measurement in the input tree
   // Make a simple tree for fitting, and then
   // call the fitter.
   cout << "HEre in MC truth" << endl;

   // The Pt and Eta histograms
   TH1F ptPass("hptpass","Pt Pass",ptNbins_,ptLow_,ptHigh_);
   TH1F ptAll("hptall","Pt All",ptNbins_,ptLow_,ptHigh_);

   TH1F etaPass("hetapass","Eta Pass",etaNbins_,etaLow_,etaHigh_);
   TH1F etaAll("hetaall","Eta All",etaNbins_,etaLow_,etaHigh_);

   //for( int i=0; i<fChain_->GetEntries(); ++i )
   for( int i=0; i<NumEvents_[0]; ++i )
   {
      fChain_->GetEntry(i);

      for( int n=0; n<ncnd_; ++n )
      {
	 if( cnd_type_[n] != tagProbeType_ ) continue;
	 
	 // These are swapped for now because of an old bug
	 if( cnd_aprobe_[n] == 1 && cnd_pprobe_[n] == 1 )
	 {
	    ptPass.Fill(cnd_pt_[n]);
	    etaPass.Fill(cnd_eta_[n]);
	 }
	 if( cnd_aprobe_[n] == 1 )
	 {
	    ptAll.Fill(cnd_pt_[n]);
	    etaAll.Fill(cnd_eta_[n]);
	 }
      }
   }

   TFile outRootFile(fitFileName_.c_str(),"RECREATE");
   outRootFile.cd();
   cout << "Writing MC Truth Eff hists!" << endl; 

   string hname = "truth_eff_Pt";
   string htitle = "Efficiency vs Pt";
   TH1F pteffhist(hname.c_str(),htitle.c_str(),ptNbins_,ptLow_,ptHigh_);
   pteffhist.Sumw2();
   pteffhist.Divide(&ptPass,&ptAll,1.0,1.0,"B");
   pteffhist.Write();

   outRootFile.cd();
   hname = "truth_eff_Eta";
   htitle = "Efficiency vs Eta";
   TH1F etaeffhist(hname.c_str(),htitle.c_str(),etaNbins_,etaLow_,etaHigh_);
   etaeffhist.Sumw2();
   etaeffhist.Divide(&etaPass,&etaAll,1.0,1.0,"B");
   etaeffhist.Write();

   outRootFile.Close();

   return;
}
// ******************************************************** //


// ------------ method called once each job just before starting event loop  ------------
void 
Histogrammer::beginJob(const edm::EventSetup&)
{

}

// ------------ method called once each job just after ending the event loop  ------------
void 
Histogrammer::endJob() {
}

//define this as a plug-in


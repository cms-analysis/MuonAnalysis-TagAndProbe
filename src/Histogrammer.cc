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
// $Id: Histogrammer.cc,v 1.3 2008/04/22 19:57:11 ahunt Exp $
//
//

#include "MuonAnalysis/TagAndProbe/interface/Histogrammer.h"
#include "MuonAnalysis/TagAndProbe/interface/RooCMSShapePdf.h"

#include <TArrow.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TGaxis.h>
#include <TGraphAsymmErrors.h>
#include <TIterator.h>
#include <TLatex.h>
#include <TString.h>
#include <TStyle.h>

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

   fileNames       = iConfig.getUntrackedParameter< vector<string> >("inputFileNames");
   quantities      = iConfig.getUntrackedParameter< vector<string> >("quantities"); 
   conditions      = iConfig.getUntrackedParameter< vector<string> >("conditions"); 
   outputFileNames = iConfig.getUntrackedParameter< vector<string> >("outputFileNames");
   XBins           = iConfig.getUntrackedParameter< vector<unsigned int> >("XBins");
   XMin            = iConfig.getUntrackedParameter< vector<double> >("XMin");
   XMax            = iConfig.getUntrackedParameter< vector<double> >("XMax");
   logY            = iConfig.getUntrackedParameter< vector<unsigned int> >("logY");
   
   // Normalization variables
   lumi_           = iConfig.getUntrackedParameter< vector<double> >("Luminosity");
   xsection_       = iConfig.getUntrackedParameter< vector<double> >("CrossSection");

   // Efficiency input variables
   calcEffsSB_     = iConfig.getUntrackedParameter< bool >("CalculateEffSideBand",false);
   calcEffsFitter_ = iConfig.getUntrackedParameter< bool >("CalculateEffFitter",false);

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
   
   bifurGaussFrac_  = iConfig.getUntrackedParameter< double >("BifurGaussFrac",0.71);

   bkgAlpha_        = iConfig.getUntrackedParameter< double >("BkgAlpha",0.67);
   bkgBeta_         = iConfig.getUntrackedParameter< double >("BkgBeta",0.05);
   bkgPeak_         = iConfig.getUntrackedParameter< double >("BkgPeak",91.1876);
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

   vectorSize = quantities.size();

   Histograms = new TH1F[vectorSize];
   NumEvents = new int[vectorSize];

   // Chain files together
   std::string tempString;
   fChain = new TChain("evttree");

   vector<string>::iterator it=fileNames.begin();
   for( int i = 0; it < fileNames.end(); it++, i++)
   {
      tempString = *it;
      fChain->Add(tempString.c_str());
      NumEvents[i] = fChain->GetEntries();
   }

   doAnalyze_ = true;
   if(outputFileNames.size() != vectorSize){
      doAnalyze_ = false;
      cout << "outputFileNames is not the same size as quantities" << endl;    
   }else if(conditions.size() != vectorSize){
      doAnalyze_ = false;
      cout << "conditions is not the same size as quantities" << endl;    
   }else if(XBins.size() != vectorSize){
      doAnalyze_ = false;
      cout << "XBins is not the same size as quantities" << endl;    
   }else if(XMax.size() != vectorSize){
      doAnalyze_ = false;
      cout << "XMax is not the same size as quantities" << endl;    
   }else if(XMin.size() != vectorSize){
      doAnalyze_ = false;
      cout << "XMin is not the same size as quantities" << endl;    
   }else if(logY.size() != vectorSize){
      doAnalyze_ = false;
      cout << "logY is not the same size as quantities" << endl;    
   }

   // Set default lumi values to 10 pb^-1
   if( lumi_.size() != fileNames.size() )
   {
      if( lumi_.size() > fileNames.size() ) lumi_.clear();
      for( int i=lumi_.size(); i<fileNames.size(); ++i ) lumi_.push_back(10.0);
   }
   // Set default cross-section values to 1000 pb
   if( xsection_.size() != fileNames.size() )
   {
      if( xsection_.size() > fileNames.size() ) xsection_.clear();
      for( int i=xsection_.size(); i<fileNames.size(); ++i ) xsection_.push_back(1000.0);
   }

   // Set the raw weights ... i.e. the number of events to
   // scale to ...
   for( int i=0; i<lumi_.size(); ++i )
   {
      if( lumi_[i]>0.0 && xsection_[i]>0.0 ) weight_.push_back(lumi_[i]*xsection_[i]);
      else                                   weight_.push_back(1.0);
   }

}

Histogrammer::~Histogrammer()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();

   delete [] Histograms;
   delete [] NumEvents;
}

// ------------ method called to for each event  ------------
void
Histogrammer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  
  
   if( doAnalyze_ )
   {
      for(unsigned int i = 0; i < vectorSize; i++)
      {
	 CreateHistogram(Histograms[i], i);
	 SaveHistogram(Histograms[i], outputFileNames[i], logY[i]);
      }
   }

   if( calcEffsSB_ || calcEffsFitter_ ) CalculateEfficiencies();
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
   tempsstream << quantities[i] << " >>  " << HistoName.str();  

   if(fChain->FindBranch(quantities[i].substr(0,quantities[i].find("[")).c_str())){ 
      Histo =  TH1F(HistoName.str().c_str(), "", XBins[i], XMin[i], XMax[i]);
      fChain->Draw(tempsstream.str().c_str());
   }else{
      std::cout << "Branch does not exist: " << quantities[i] << std::endl;
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

   const Int_t PeakBin = (Int_t)(Peak - xmin)/BinWidth + 1; // Peak
   const Double_t SDBin = SD/BinWidth; // Standard deviation
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

   TH1F TPMassHistoSig("TPMassHistoSig", "", 200,0,200);
   TH1F TPMassHistoBG("TPMassHistoBG", "", 200,0,200);
   TH1F TPMassHistoTot("TPMassHistoTot", "", 200,0,200);
   TH1F TPMassHistoSBS("TPMassHistoSBS","", 200,0,200);

   fChain->Draw("tp_mass >> TPMassHistoSig","", "", NumEvents[0], 0);
   fChain->Draw("tp_mass >> TPMassHistoBG", "", "", NumEvents[1] - NumEvents[0], NumEvents[0] + 1);
  
   const Double_t SigScale = weight_[0]/NumEvents[0];

   std::cout << "SigScale: " << SigScale << std::endl;

   TPMassHistoTot.Add(&TPMassHistoSig, &TPMassHistoBG, SigScale, 1);
   SideBandSubtraction(TPMassHistoTot, TPMassHistoSBS, 90, 2);

   SaveHistogram(TPMassHistoSig, "TPMassHistoSig.png", 1);
   SaveHistogram(TPMassHistoBG, "TPMassHistoBG.png", 1);
   SaveHistogram(TPMassHistoTot, "TPMassHistoTot.png", 1);
   SaveHistogram(TPMassHistoSBS, "TPMassHistoSBS.png", 1);

}
// ******************************************** //


// ********** Z -> l+l- Fitter ********** //
void Histogrammer::ZllEffFitter( TTree &fitTree, string &fileName, string &bvar,
				 int bnbins, double blow, double bhigh )
{
   string fmode = "RECREATE";
   if( bvar == "Eta" ) fmode = "UPDATE";
   TFile outRootFile(fileName.c_str(),fmode.c_str());
   outRootFile.cd();
   //fitTree.Write();
   
   string hname = "heff_" + bvar;
   string htitle = "Efficiency vs " + bvar;
   TH1F effhist(hname.c_str(),htitle.c_str(),bnbins,blow,bhigh);

   double bwidth = (bhigh-blow)/(double)bnbins;

   for( int bin=0; bin<bnbins; ++bin )
   {

      // The fit variable - lepton invariant mass
      RooRealVar Mass("Mass","Invariant Lepton Mass", massLow_, massHigh_, "GeV");

      // The binning variable
      string bunits = "";
      double lowEdge = blow + (double)bin*bwidth;
      double highEdge = lowEdge + bwidth;
      if( bvar == "Pt" ) bunits = "GeV";
      RooRealVar Pt(bvar.c_str(),bvar.c_str(),lowEdge,highEdge,bunits.c_str());

      // The weighting
      RooRealVar Weight("Weight","Weight",0.0,10000.0);
 
      // ********** Construct signal shape PDF ********** //

      // Signal PDF variables
      RooRealVar signalMean("signalMean","signalMean",
			    signalMean_[0],signalMean_[1],signalMean_[2]);
      RooRealVar signalWidth("signalWidth","signalWidth",
			     signalWidth_[0],signalWidth_[1],signalWidth_[2]);
      RooRealVar signalSigma("signalSigma","signalSigma",
			     signalSigma_[0],signalSigma_[1],signalSigma_[2]);
      RooRealVar signalWidthL("signalWidthL","signalWidthL",
			      signalWidthL_[0],signalWidthL_[1],signalWidthL_[2]);
      RooRealVar signalWidthR("signalWidthR","signalWidthR",
			      signalWidthR_[0],signalWidthR_[1],signalWidthR_[2]);
  
      // Voigtian
      RooVoigtian signalVoigtPdf("signalVoigtPdf", "signalVoigtPdf", 
				 Mass, signalMean, signalWidth, signalSigma);

      // Bifurcated Gaussian
      RooBifurGauss signalGaussBifurPdf("signalGaussBifurPdf", "signalGaussBifurPdf", 
					Mass, signalMean, signalWidthL, signalWidthR);

      // Bifurcated Gaussian fraction
      RooRealVar fBifurGauss("fBifurGauss","fBifurGauss",bifurGaussFrac_);

      // The total signal PDF
      RooAddPdf  signalShapePdf("signalShapePdf", "signalShapePdf",
				signalVoigtPdf,signalGaussBifurPdf,fBifurGauss);

      // ********** Construct background shape PDF ********** //

      // Background PDF variables
      RooRealVar bkgAlpha("bkgAlpha","bkgAlpha",bkgAlpha_);
      RooRealVar bkgBeta("bkgBeta","bkgBeta",bkgBeta_);
      RooRealVar bkgGamma("bkgGamma","bkgGamma",bkgGamma_[0],bkgGamma_[1],bkgGamma_[2]);
      RooRealVar bkgPeak("bkgPeak","bkgPeak",bkgPeak_);

      // CMS Background shape
      RooCMSShapePdf bkgShapePdf("bkgShapePdf","bkgShapePdf", 
				 Mass,bkgAlpha,bkgBeta,bkgGamma,signalMean) ;

      // Now define some efficiency variables  
      RooRealVar efficiency("efficiency","efficiency",
			    efficiency_[0],efficiency_[1],efficiency_[2]);
      RooRealVar nSig("nSig","nSig",
		      numSignal_[0],numSignal_[1],numSignal_[2]);
      RooRealVar nBkgpass("nBkgpass","nBkgpass",
			  numBkgPass_[0],numBkgPass_[1],numBkgPass_[2]);
      RooRealVar nBkgfail("nBkgfail","nBkgfail",
			  numBkgFail_[0],numBkgFail_[1],numBkgFail_[2]);

      RooFormulaVar nSigpass("nSigpass","nSig*efficiency", RooArgList(nSig,efficiency) );
      RooFormulaVar nSigfail("nSigfail","nSig*(1.0 - efficiency)", RooArgList(nSig,efficiency) );

      RooArgList componentspass(signalShapePdf,bkgShapePdf);
      RooArgList yieldspass(nSigpass, nBkgpass);
      RooArgList componentsfail(signalShapePdf,bkgShapePdf);
      RooArgList yieldsfail(nSigfail, nBkgfail);	  

      RooAddPdf sumpass("sumpass","fixed extended sum pdf",componentspass,yieldspass);
      RooAddPdf sumfail("sumfail","fixed extended sum pdf",componentsfail, yieldsfail);

      // Make the category variable that defines the two fits,
      // namely whether the probe passes or fails the eff criteria.
      RooCategory ProbePass("ProbePass","sample");
      ProbePass.defineType("pass",1);
      ProbePass.defineType("fail",0);  
  
      // The total simultaneous fit ...
      RooSimultaneous totalPdf("totalPdf","totalPdf",ProbePass);
      ProbePass.setLabel("pass");
      totalPdf.addPdf(sumpass,ProbePass.getLabel());
      totalPdf.Print();
      ProbePass.setLabel("fail");
      totalPdf.addPdf(sumfail,ProbePass.getLabel());
      totalPdf.Print();

      // Add the TTree as our data set ... with the weight in case 
      // we are using chained MC
      RooDataSet* data = new RooDataSet("fitData","fitData",&fitTree,
					RooArgSet(ProbePass,Mass,Pt,Weight));
      data->get()->Print();
      data->setWeightVar("Weight");
      data->get()->Print();

      RooDataHist *bdata = new RooDataHist("bdata","Binned Data",
					   RooArgList(Mass,ProbePass),*data);

      // Count the number of passing and failing probes in the region
      // making sure we have enough to fit ...
      cout << "About to count the number of events" << endl;
      int npassR = data->sumEntries("ProbePass==1");
      int nfailR = data->sumEntries("ProbePass==0");
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
	    cout << "GOF Entries " << dset->numEntries() << " " <<type->GetName() << std::endl;
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
      //RooNLLVar nll("nll","nll",totalPdf,*data,kTRUE);
      //RooMinuit m(nll);
      RooChi2Var chi2("chi2","chi2",totalPdf,*bdata,DataError(RooAbsData::SumW2));
      RooMinuit m(chi2);
      m.setErrorLevel(0.5); // <<< HERE
      //m.setStrategy(2);
      //m.hesse();
      m.migrad();
      //m.hesse();
      //m.minos();
      fitResult = m.save();

      std::cout << "Signal yield: " << nSig.getVal() << " +- "
		<< nSig.getError() << " + " << nSig.getAsymErrorHi()
		<<" - "<< nSig.getAsymErrorLo() << std::endl;
      std::cout << "Efficiency: "<< efficiency.getVal() << " +- "
		<< efficiency.getError() << " + " << efficiency.getAsymErrorHi()
		<<" + "<< efficiency.getAsymErrorLo() << std::endl;

      // Fill the efficiency hist
      effhist.SetBinContent(bin+1,efficiency.getVal());
      effhist.SetBinError(bin+1,efficiency.getError());

      // ********** Make and save Canvas for the plots ********** //
      outRootFile.cd();
      
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

      RooPlot* frame1 = Mass.frame(60);
      frame1->SetTitle("");
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
      RooPlot* frame2 = Mass.frame(60);
      frame2->SetTitle("");
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
      RooPlot* frame3 = Mass.frame(60);
      frame3->SetTitle("");
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
   if( calcEffsFitter_ )
   {
      // Loop over the number of different types of 
      // efficiency measurement in the input tree
      // Make a simple tree for fitting, and then
      // call the fitter.
      int tagprobe_type = 0;

      Int_t           nrtp;
      //Int_t           tp_type[100];
      Int_t           tp_true[100];
      Int_t           tp_ppass[100];
      Float_t         tp_mass[100]; 
      Float_t         tp_dpt[100][2];
      Float_t         tp_deta[100][2];
      
      TBranch        *b_nrtp;   
      TBranch        *b_tp_true;   
      //TBranch        *b_tp_type;   
      TBranch        *b_tp_ppass;   
      TBranch        *b_tp_mass;   
      TBranch        *b_tp_dpt;   
      TBranch        *b_tp_deta;   

      fChain->SetBranchAddress("nrtp", &nrtp, &b_nrtp);
      fChain->SetBranchAddress("tp_true", tp_true, &b_tp_true);
      //fChain->SetBranchAddress("tp_type", tp_type, &b_tp_type);
      fChain->SetBranchAddress("tp_ppass", tp_ppass, &b_tp_ppass);
      fChain->SetBranchAddress("tp_mass", tp_mass, &b_tp_mass);
      fChain->SetBranchAddress("tp_dpt", tp_dpt, &b_tp_dpt);
      fChain->SetBranchAddress("tp_deta", tp_deta, &b_tp_deta);

      fChain->SetBranchStatus("*",0);
      fChain->SetBranchStatus("nrtp",1);
      //fChain->SetBranchStatus("tp_type",1);
      fChain->SetBranchStatus("tp_ppass",1);
      fChain->SetBranchStatus("tp_mass",1);
      fChain->SetBranchStatus("tp_dpt",1);
      fChain->SetBranchStatus("tp_deta",1);

      // Make the simple fit tree
      int    ProbePass;
      double Mass;
      double Pt;
      double Eta;
      double Weight;

      TTree fitTree("fitter_tree","Tree For Fitting",1);
      fitTree.Branch("ProbePass",&ProbePass,"ProbePass/I");
      fitTree.Branch("Mass",     &Mass,     "Mass/D");
      fitTree.Branch("Pt",       &Pt,       "Pt/D");
      fitTree.Branch("Eta",      &Eta,      "Eta/D");
      fitTree.Branch("Weight",   &Weight,   "Weight/D");

      int nFile = 0;
      Weight = 1.0;
      if( NumEvents[nFile] > 0 ) Weight = weight_[nFile]/(double)NumEvents[nFile];
      cout << "Filling fit tree with weight " << Weight << endl;
      for( int i=0; i<fChain->GetEntries(); ++i )
      {
	 fChain->GetEntry(i);

	 // If we pass a boundary, change the weight
	 if( nFile<(int)fileNames.size() && i>=NumEvents[nFile] )
	 {
	    ++nFile;
	    if( lumi_[nFile] > 0.0 )
	    {
	       if( NumEvents[nFile] > 0 ) Weight = weight_[nFile]/(double)NumEvents[nFile];
	    }
	    else
	    {
	       Weight = weight_[nFile];
	    }
	    cout << "Filling fit tree with weight " << Weight << endl;
	 }

	 for( int n=0; n<nrtp; ++n )
	 {
	    //if( tp_type[n] != tagprobe_type ) continue;

	    ProbePass = tp_ppass[n];
	    Mass      = (double)tp_mass[n];
	    Pt        = (double)tp_dpt[n][1];
	    Eta       = (double)tp_deta[n][1];

	    fitTree.Fill();
	 }
      }

      // We have filled the simple tree ... call the fitter
      string binnedVar = "Pt";
      ZllEffFitter( fitTree, fitFileName_, binnedVar, ptNbins_, ptLow_, ptHigh_ );
      binnedVar = "Eta";
      ZllEffFitter( fitTree, fitFileName_, binnedVar, etaNbins_, etaLow_, etaHigh_ );

   }

   if( calcEffsSB_ )
   {
      // Need to fill this function
      return;
   }
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


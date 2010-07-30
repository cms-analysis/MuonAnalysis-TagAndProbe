#include <memory>
#include <string>
#include <fstream>
#include <iostream>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CondCore/DBOutputService/interface/PoolDBOutputService.h"
#include "CondFormats/PhysicsToolsObjects/interface/PerformancePayloadFromTable.h"

#include "CondFormats/PhysicsToolsObjects/interface/PerformanceWorkingPoint.h"

#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TROOT.h>
#include <RooDataSet.h>
#include <RooRealVar.h>
#include <RooDataHist.h>
#include <RooArgSet.h>

using namespace std;
using namespace RooFit;

class PhysicsPerformanceDBWriterFromTPDataset : public edm::EDAnalyzer
{
public:
  PhysicsPerformanceDBWriterFromTPDataset(const edm::ParameterSet&);
  virtual void beginJob();
  virtual void analyze(const edm::Event&, const edm::EventSetup&) {}
  virtual void endJob() {}
  ~PhysicsPerformanceDBWriterFromTPDataset() {}

private:
  std::vector<std::string> rec1;
  std::vector<std::string> rec2;
  boost::uint64_t iovBegin;
  boost::uint64_t iovEnd;  
  std::vector<std::string> inputHistoFiles;
  std::vector<std::string> inputDatasetNames; 
  std::vector<std::string> inputBinningVariables;
  std::vector<std::string> inputResultTypes;
  std::vector<std::string> inputAlgorithmNames;
  std::vector<double> inputDiscriminatorCuts;
  std::string inputConcreteClass; 

};

PhysicsPerformanceDBWriterFromTPDataset::PhysicsPerformanceDBWriterFromTPDataset
  (const edm::ParameterSet& p)
{
  rec1 = p.getParameter< std::vector <std::string> >("RecordPayloads");
  rec2 = p.getParameter< std::vector <std::string> >("RecordWPs");
  iovBegin = p.getParameter<boost::uint64_t>("IOVBegin"); 
  iovEnd = p.getParameter<boost::uint64_t>("IOVEnd"); 
  inputHistoFiles = p.getParameter< std::vector<std::string> >("inputHistoFiles");
  inputDatasetNames = p.getParameter< std::vector<std::string> >("inputDatasetNames"); 
  inputAlgorithmNames = p.getParameter< std::vector<std::string> >("inputAlgorithmNames");
  inputDiscriminatorCuts = p.getParameter< std::vector<double> >("inputDiscriminatorCuts");
  inputConcreteClass = p.getUntrackedParameter<std::string>("inputConcreteClass","PerformancePayloadFromTable");
  inputResultTypes = p.getParameter< std::vector<std::string> >("inputResultTypes");
  inputBinningVariables = p.getParameter< std::vector<std::string> >("inputBinningVariables");  

}

void PhysicsPerformanceDBWriterFromTPDataset::beginJob()
{
  //
  // read object from file
  //

  //
  // File Format is
  // - tagger name
  // - cut
  // - concrete class name
  // - how many results and how many binning
  // - their values
  // - vector<float>
  //

  std::string tagger; 
  float cut;  
  std::string concreteType; 
  std::string comment; 
  std::vector<float> pl; 
  std::string infilename;
  std::string datasetname;
  int stride; 
  int nres, nbin;  
  int number = 0;
  std::string tmprec1, tmprec2;
  boost::uint64_t tmpiovBegin, tmpiovEnd;

  std::vector<PerformanceResult::ResultType> res;   
  std::vector<BinningVariables::BinningVariablesType> bin; 

  concreteType = inputConcreteClass;  

  nres = inputResultTypes.size();
  cout << nres << " result types" << endl;

  // First read the result types and binning variables
  //  for(unsigned int i = 0; i < inputResultTypes.size(); ++i)
  //    {
      //      int tmp = -1; 
      //      if(inputResultTypes[i] == "efficiency")
      //	tmp = 2001;
  //      else if(inputResultTypes[i] == "efficiency_symerr")
  //	tmp = 2002; 
  //      else
  //	cout << "Unknown efficeincy result type " << inputResultTypes[i] << endl;
  //      cout << "\tResult type = " << tmp << " (" << inputResultTypes[i] << ")" << endl;
  //      res.push_back((PerformanceResult::ResultType)(tmp)); 
  //    } 

  int tmp = 2001;
  res.push_back((PerformanceResult::ResultType)(tmp));
  tmp = 2002;
  res.push_back((PerformanceResult::ResultType)(tmp));

  nbin = inputBinningVariables.size();
  cout << nbin << " binning variables" << endl;

  for(unsigned int i = 0; i < inputBinningVariables.size(); ++i) 
    {
      int tmp = -1; 
      if(inputBinningVariables[i] == "pt")
	tmp = 1001;
      else if(inputBinningVariables[i] == "charge")
	tmp = 1002;
      else if(inputBinningVariables[i] == "abseta" || inputBinningVariables[i] == "eta")
	tmp = 1003;
      else if(inputBinningVariables[i] == "phi")
	tmp = 1004;
      else 
	cout << "Unknown binning variable type " << inputBinningVariables[i] << endl; 
      cout << "\tBinning variable = " << tmp << " (" << inputBinningVariables[i] << ")" << endl; 

      bin.push_back((BinningVariables::BinningVariablesType)(tmp)); 
    } 

  stride = nres+nbin*2; 

  // Now read the acutal payload
  for(int m = 0; m < 3; m++)
    {
  for(unsigned int i = 0; i < inputHistoFiles.size(); ++i)
    {
      number = 0; 
      pl.clear();

      infilename = inputHistoFiles[i]; 
      datasetname = inputDatasetNames[i]; 
      tagger = inputAlgorithmNames[i];
      cut = inputDiscriminatorCuts[i];    
      tmpiovBegin = iovBegin;
      tmpiovEnd = iovEnd;
      tmprec1 = rec1[i];
      tmprec2 = rec2[i];

      if(m == 1)
	{
	  tmprec1 = "LOERR_" + tmprec1;
          tmprec2 = "LOERR_" + tmprec2; 
	  tagger = tagger + "_LowerError";
	}
      if(m == 2)
	{
	  tmprec1 = "UPERR_" + tmprec1;
          tmprec2 = "UPERR_" + tmprec2; 
	  tagger = tagger + "_UpperError";
	}

      cout << "Reading from Tag-&-Probe file " << infilename << endl; 
      cout << "\tReading efficiencies from RooDataset named " << datasetname << endl;
      cout << "Algorithm name = " << tagger << endl;
      cout << "Discriminator cut = " << cut << endl;

      RooDataSet *datatmp = 0;

      TFile *f = TFile::Open(infilename.c_str());
      f->cd();

      datatmp = (RooDataSet *)f->Get(datasetname.c_str());
      const RooArgSet* vars = datatmp->get();
      RooRealVar* effvar = (RooRealVar*) vars->find("efficiency");
      RooRealVar* etavar = (RooRealVar*) vars->find("abseta");
      RooRealVar* ptvar = (RooRealVar*) vars->find("pt");

      // Fill DB tables
      float binlowedgex = -999.;
      float binhighedgex = -999.;
      float binlowedgey = -999.; 
      float binhighedgey = -999.; 
      float bincenterx = -999.;
      float bincentery = -999.;
      float bincontent = -999.;
      float binerror = -999.;

      TH2F *h = new TH2F("eff", "eff", ptvar->getBinning().numBins(),
			 ptvar->getBinning().array(), etavar->getBinning().numBins(),
			 etavar->getBinning().array());
      TH2F* hlo = (TH2F *)h->Clone("eff_lo");
      TH2F* hhi = (TH2F *)h->Clone("eff_hi");

      int l = 0;

      for(int j = 0;j < ptvar->getBinning().numBins();j++)
	{
	  for(int k = 0;k < etavar->getBinning().numBins();k++)
	    {
	      binlowedgex = ptvar->getBinning().binLow(j);
	      binhighedgex = ptvar->getBinning().binHigh(j);
	      binlowedgey = etavar->getBinning().binLow(k);
	      binhighedgey = etavar->getBinning().binHigh(k) ;
	      bincenterx = ptvar->getBinning().binCenter(j);
	      bincentery = etavar->getBinning().binCenter(k);
	      int b = h->FindBin(bincenterx, bincentery);

	      datatmp->get(l);
	      h->SetBinContent(b, effvar->getVal());
	      hlo->SetBinContent(b, effvar->getVal()+effvar->getErrorLo());
	      hhi->SetBinContent(b, effvar->getVal()+effvar->getErrorHi());

	      if(m == 0)
		bincontent = effvar->getVal();
	      if(m == 1)
		bincontent = effvar->getErrorLo();
	      if(m == 2)
		bincontent = effvar->getErrorHi();

	      //	      binerror = fabs(effvar->getErrorHi() - effvar->getErrorLo())/2.0;
	      binerror = 0.0;

              pl.push_back(binlowedgex);
              pl.push_back(binhighedgex);
              pl.push_back(binlowedgey);
              pl.push_back(binhighedgey);

	      cout << "(" << binlowedgex << " < pT < "  << binhighedgex << "), "
		   << "(" << binlowedgey << " < eta < " << binhighedgey << ") "
		   << " = " << h->GetBinContent(b)
		   << "+" << (hhi->GetBinContent(b) - h->GetBinContent(b))
		   << "-" << (h->GetBinContent(b) - hlo->GetBinContent(b)) << endl;

              cout << " Inserting " << bincontent << " in position " << number
                   << " (" << binlowedgex << " < " << inputBinningVariables[0] << " < " << binhighedgex << ", "
                   << binlowedgey << " < " << inputBinningVariables[1] << " < " << binhighedgey << ")" << endl;
              pl.push_back(bincontent);
              number++;
              cout << " Inserting " << binerror << " in position " << number
                   << " (" << binlowedgex << " < " << inputBinningVariables[0] << " < " << binhighedgex << ", "
                   << binlowedgey << " < " << inputBinningVariables[1] << " < " << binhighedgey << ")" << endl;
              pl.push_back(binerror);

              number++;
	      l++;
	    }
	}

      f->Close();

      if (stride != nbin*2+nres)
	{ 
	  std::cout <<" Table not well formed"<<std::endl; 
	} 
      if ((number % stride) != 0)
	{ 
	  std::cout <<" Table not well formed"<<std::endl; 
	} 

      PerformanceWorkingPoint * wp = new PerformanceWorkingPoint(cut, tagger); 
      
      PerformancePayloadFromTable * mupl = 0; 
      
      if (concreteType == "PerformancePayloadFromTable")
	{ 
	  mupl = new PerformancePayloadFromTable(res, bin, stride, pl); 
	}
      else
	{ 
	  cout <<" Non existing request: " <<concreteType<<endl; 
	} 
   
      cout <<" Created the "<<concreteType <<" object"<<endl; 

      // 
      // now create pl etc etc 
      // 
 
      edm::Service<cond::service::PoolDBOutputService> s; 
      if (s.isAvailable()) 
	{ 
	  if (s->isNewTagRequest(tmprec1)) 
	    { 
	      s->createNewIOV<PerformancePayload>(mupl, 
						  tmpiovBegin, 
						  tmpiovEnd, 
						  tmprec1); 
	    } 
	  else 
	    { 
	      s->appendSinceTime<PerformancePayload>(mupl, 
						     tmpiovBegin, 
						     tmprec1); 
	    } 
	  cout << "Wrote payload" << endl;
	} 

      // write also the WP 
   
      if (s.isAvailable()) 
	{ 
	  if (s->isNewTagRequest(tmprec2)) 
	    { 
	      s->createNewIOV<PerformanceWorkingPoint>(wp, 
						       tmpiovBegin, 
						       tmpiovEnd, 
						       tmprec2); 
	    } 
	  else 
	    { 
           
	      s->appendSinceTime<PerformanceWorkingPoint>(wp, 
							  tmpiovBegin, 
							  tmprec2); 
	    } 
	  cout << "Wrote Working Point" << endl;
	}       

    }
    }
}

DEFINE_FWK_MODULE(PhysicsPerformanceDBWriterFromTPDataset);

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
  unsigned long long iovBegin;
  unsigned long long iovEnd;  
  std::vector<std::string> inputHistoFiles;
  std::vector<std::string> inputDatasetNames; 
  std::vector<std::string> inputBinningVariables;
  std::vector<std::string> inputResultTypes;
  std::vector<std::string> inputAlgorithmNames;
  std::vector<double> inputDiscriminatorCuts;
  std::string inputConcreteClass; 
  bool inputMergeTwoInputDatasets;
  double inputMergeDatasetsPtBoundary; 

};

PhysicsPerformanceDBWriterFromTPDataset::PhysicsPerformanceDBWriterFromTPDataset
  (const edm::ParameterSet& p)
{
  rec1 = p.getParameter< std::vector <std::string> >("RecordPayloads");
  rec2 = p.getParameter< std::vector <std::string> >("RecordWPs");
  iovBegin = p.getParameter<unsigned long long>("IOVBegin"); 
  iovEnd = p.getParameter<unsigned long long>("IOVEnd"); 
  inputHistoFiles = p.getParameter< std::vector<std::string> >("inputHistoFiles");
  inputDatasetNames = p.getParameter< std::vector<std::string> >("inputDatasetNames"); 
  inputAlgorithmNames = p.getParameter< std::vector<std::string> >("inputAlgorithmNames");
  inputDiscriminatorCuts = p.getParameter< std::vector<double> >("inputDiscriminatorCuts");
  inputConcreteClass = p.getUntrackedParameter<std::string>("inputConcreteClass","PerformancePayloadFromTable");
  inputResultTypes = p.getParameter< std::vector<std::string> >("inputResultTypes");
  inputBinningVariables = p.getParameter< std::vector<std::string> >("inputBinningVariables");  
  inputMergeTwoInputDatasets = p.getParameter< bool >("inputMergeTwoInputDatasets"); 
  inputMergeDatasetsPtBoundary = p.getParameter< double >("inputMergeDatasetsPtBoundary");
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
  std::string infilename, infilename2;
  std::string datasetname, datasetname2;
  int stride; 
  int nres, nbin;  
  int number = 0;
  std::string tmprec1, tmprec2;
  unsigned long long tmpiovBegin, tmpiovEnd;

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

      cout << "Reading from Tag-&-Probe file " << infilename << endl;
      cout << "\tReading efficiencies from RooDataset named " << datasetname << endl;

      // Allow for reading from two datasets, and loading into a single table 
      // (e.g. for JPsi+Z based efficiencies) 
      if(inputMergeTwoInputDatasets == true)
	{
	  infilename2 = inputHistoFiles[i+1];
	  datasetname2 = inputDatasetNames[i+1];
	  cout << "Reading from Tag-&-Probe file " << infilename2 << endl;
	  cout << "\tReading efficiencies from RooDataset named " << datasetname2 << endl;
	}

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

      cout << "\tPayload " << tmprec1 << endl; 
      cout << "\tWP " << tmprec2 << endl;
      cout << "\tAlgorithm name = " << tagger << endl;
      cout << "\tDiscriminator cut = " << cut << endl;

      RooDataSet *datatmp = 0;

      TFile *f = 0;
      f = TFile::Open(infilename.c_str());
      f->cd();

      datatmp = (RooDataSet *)f->Get(datasetname.c_str());
      const RooArgSet* vars = datatmp->get();
      RooRealVar* effvar = (RooRealVar*) vars->find("efficiency");
      RooRealVar* etavar = (RooRealVar*) vars->find("abseta");
      RooRealVar* ptvar = (RooRealVar*) vars->find("pt");

      RooDataSet *datatmp2 = 0;
      TFile *f2 = 0;
      RooRealVar* effvar2 = 0;
      RooRealVar* etavar2 = 0;
      RooRealVar* ptvar2 = 0;

      if(inputMergeTwoInputDatasets == true)
	{
	  f2 = TFile::Open(infilename2.c_str());

	  datatmp2 = (RooDataSet *)f2->Get(datasetname2.c_str());
	  const RooArgSet* vars2 = datatmp2->get();
	  effvar2 = (RooRealVar*) vars2->find("efficiency");
	  etavar2 = (RooRealVar*) vars2->find("abseta");
	  ptvar2 = (RooRealVar*) vars2->find("pt");
	}

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

      TH2F *h2 = 0;
      TH2F* hlo2 = 0;
      TH2F* hhi2 = 0;

      if(inputMergeTwoInputDatasets == true)
	{
	  h2 = new TH2F("eff2", "eff2", ptvar2->getBinning().numBins(),
			     ptvar2->getBinning().array(), etavar2->getBinning().numBins(),
			     etavar2->getBinning().array());
	  cout << "Cloning for efficiency bin boundaries" << endl;
	  hlo2 = (TH2F *)h2->Clone("eff2_lo");
	  hhi2 = (TH2F *)h2->Clone("eff2_hi");
	}

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

	      if((inputMergeTwoInputDatasets == false) || 
		 ((inputMergeTwoInputDatasets == true) && (binhighedgex <= inputMergeDatasetsPtBoundary)))
		{
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
	}

      f->Close();

      if(inputMergeTwoInputDatasets == true)
	{
	  f2->cd();
	  int l = 0;

	  for(int j = 0;j < ptvar2->getBinning().numBins();j++)
	    {
	      for(int k = 0;k < etavar2->getBinning().numBins();k++)
		{
		  binlowedgex = ptvar2->getBinning().binLow(j);
		  binhighedgex = ptvar2->getBinning().binHigh(j);
		  binlowedgey = etavar2->getBinning().binLow(k);
		  binhighedgey = etavar2->getBinning().binHigh(k) ;
		  bincenterx = ptvar2->getBinning().binCenter(j);
		  bincentery = etavar2->getBinning().binCenter(k);
		  int b = h2->FindBin(bincenterx, bincentery);

		  datatmp2->get(l);
		  h2->SetBinContent(b, effvar2->getVal());
		  hlo2->SetBinContent(b, effvar2->getVal()+effvar2->getErrorLo());
		  hhi2->SetBinContent(b, effvar2->getVal()+effvar2->getErrorHi());

		  if(m == 0)
		    bincontent = effvar2->getVal();
		  if(m == 1)
		    bincontent = effvar2->getErrorLo();
		  if(m == 2)
		    bincontent = effvar2->getErrorHi();

		  //              binerror = fabs(effvar2->getErrorHi() - effvar2->getErrorLo())/2.0;
		  binerror = 0.0;

		  if(binlowedgex >= inputMergeDatasetsPtBoundary)
		    {
		      pl.push_back(binlowedgex);
		      pl.push_back(binhighedgex);
		      pl.push_back(binlowedgey);
		      pl.push_back(binhighedgey);
		      
		      cout << "(" << binlowedgex << " < pT < "  << binhighedgex << "), "
			   << "(" << binlowedgey << " < eta < " << binhighedgey << ") "
			   << " = " << h2->GetBinContent(b)
			   << "+" << (hhi2->GetBinContent(b) - h2->GetBinContent(b))
			   << "-" << (h2->GetBinContent(b) - hlo2->GetBinContent(b)) << endl;
		      
		      cout << " Inserting " << bincontent << " in position " << number
			   << " (" << binlowedgex << " < " << inputBinningVariables[0] << " < " << binhighedgex << ", "
			   << binlowedgey << " < " << inputBinningVariables[1] << " < " << binhighedgey << ")" << endl;
		      
		      pl.push_back(bincontent);

		      cout << " Inserting " << binerror << " in position " << number 
			   << " (" << binlowedgex << " < " << inputBinningVariables[0] << " < " << binhighedgex << ", " 
			   << binlowedgey << " < " << inputBinningVariables[1] << " < " << binhighedgey << ")" << endl; 
		      pl.push_back(binerror); 

		      number++;
		      l++;
		    }
		}
	    }
	  i++;
	  f2->Close();
	}	  

      if (stride != nbin*2+nres)
	{ 
	  std::cout <<" \tTable not well formed"<<std::endl; 
	} 
      if ((number % stride) != 0)
	{ 
	  std::cout <<"\tTable not well formed"<<std::endl; 
	} 

      PerformanceWorkingPoint * wp = new PerformanceWorkingPoint(cut, tagger); 
      
      PerformancePayloadFromTable * mupl = 0; 
      
      if (concreteType == "PerformancePayloadFromTable")
	{ 
	  mupl = new PerformancePayloadFromTable(res, bin, stride, pl); 
	}
      else
	{ 
	  cout <<"\tNon existing request: " <<concreteType<<endl; 
	} 
   
      cout <<"\tCreated the "<<concreteType <<" object"<<endl; 

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
	  cout << "\tWrote payload" << endl;
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
	  cout << "\tWrote Working Point" << endl;
	}       

    }
    }
}

DEFINE_FWK_MODULE(PhysicsPerformanceDBWriterFromTPDataset);

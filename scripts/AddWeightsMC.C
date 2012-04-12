void AddWeightsMC(){

  //Run2011A Reweighting -- Muons
  //TFile *inputData = TFile::Open("tnpZ_PD_SingleMu_Run2011A_Total.root");
  //TFile *inputMC   = TFile::Open("tnpZ_MC_Merge_MuMu_2.1invfb_AllFall11_2011AReweighted.root","update");

  //Run2011B Reweighting -- Muons
  //TFile *inputData = TFile::Open("tnpZ_PD_SingleMu_Run2011B_Total.root");
  //TFile *inputMC   = TFile::Open("tnpZ_MC_Merge_MuMu_2.5invfb_AllFall11_2011BReweighted.root","update");

  //Run2011AB Reweighting -- Muons
  //TFile *inputData = TFile::Open("tnpZ_PD_SingleMu_Run2011AB_TOTAL.root");
  //TFile *inputMC   = TFile::Open("tnpZ_MC_Merge_MuMu_2.5invfb_AllFall11_All2011Reweighted.root","update");

  //Run2011A Reweighting -- Electrons
  //TFile *inputData = TFile::Open("tnpZ_PD_DoubleElectron_Run2011A_Total.root");
  //TFile *inputMC   = TFile::Open("tnpZ_MC_Merge_EE_2.1invfb_AllFall11_2011AReweighted.root","update");

  //Run2011B Reweighting -- Electrons
  //TFile *inputData = TFile::Open("tnpZ_PD_DoubleElectron_Run2011B_Total.root");
  //TFile *inputMC   = TFile::Open("tnpZ_MC_Merge_EE_2.5invfb_AllFall11_2011BReweighted.root","update");

  //Run2011AB Reweighting -- Electrons
  TFile *inputData = TFile::Open("tnpZ_PD_DoubleElectron_Run2011AB_TOTAL.root");
  TFile *inputMC   = TFile::Open("tnpZ_MC_Merge_EE_2.5invfb_AllFall11_All2011Reweighted.root","update");

  Float_t DatanVtx;
  Float_t MCnVtx;

  Float_t weight;

  Float_t MCMass;
  Float_t MCabseta;
  Float_t MCeta;
  Int_t MC_passTM;
  Int_t MC_passID;
  Int_t MC_passIP;
  Int_t MC_passIso;

  TH1F* hDatanVtx_aux   = new TH1F("hDatanVtx_aux"  ,"",50,0.,50.);
  TH1F* hDatanVtx2_aux  = new TH1F("hDatanVtx2_aux" ,"",50,0.,50.);
  TH1F* hMCnVtx_aux     = new TH1F("hMCnVtx_aux"  ,"",50,0.,50.);
  TH1F* hWeightnVtx_aux = new TH1F("hWeightnVtx_aux","",50,0.,50.);

  TH1F* hDatanVtx   = new TH1F("hDatanVtx"  ,"",50,0.,50.);
  TH1F* hDatanVtx2  = new TH1F("hDatanVtx2" ,"",50,0.,50.);

  TH1F* hMCnVtx     = new TH1F("hMCnVtx"    ,"",50,0.,50.);
  TH1F* hWeightnVtx = new TH1F("hWeightnVtx","",50,0.,50.);

  TH1F* hWeightedMC = new TH1F("hWeightedMC","",50,0.,50.);

  hDatanVtx->Sumw2();
  hMCnVtx->Sumw2();

  hDatanVtx_aux->Sumw2();
  hDatanVtx2_aux->Sumw2();
  hMCnVtx_aux->Sumw2();

  hWeightedMC->Sumw2();

  //inputData->cd("tpTreeMuMu");
  inputData->cd("tpTreeElElBDT");
  fitter_tree->SetBranchAddress("nVtx",&DatanVtx);
  Int_t nData = fitter_tree->GetEntries();

  Double_t normData = 1.0/nData;

  for(size_t nEvent=0; nEvent<nData; ++nEvent){

    fitter_tree->GetEntry(nEvent);
    
    hDatanVtx->Fill(DatanVtx);
    hDatanVtx2->Fill(DatanVtx);

    hDatanVtx_aux->Fill(DatanVtx,normData);
    hDatanVtx2_aux->Fill(DatanVtx,normData);

  }

  //inputMC->cd("tpTreeMuMu");
  inputMC->cd("tpTreeElElBDT");
  fitter_tree->SetBranchAddress("nVtx",&MCnVtx);
  fitter_tree->SetBranchAddress("mass",&MCMass);
  fitter_tree->SetBranchAddress("abseta",&MCabseta);
  fitter_tree->SetBranchAddress("eta",&MCeta);
  fitter_tree->SetBranchAddress("passTM",&MC_passTM);
  fitter_tree->SetBranchAddress("passIP",&MC_passIP);
  fitter_tree->SetBranchAddress("passID",&MC_passID);
  fitter_tree->SetBranchAddress("passIso",&MC_passIso);
  Int_t nMC = fitter_tree->GetEntries();


  //Double_t ratio = ( (double)nData)/nMC;

  //std::cout<<" nData/nMC = "<< nData << "/" << nMC << std::endl;
  //std::cout<<" nData/nMC = "<< ratio << std::endl;

  Double_t normMC = 1.0/nMC;

  for(size_t nEvent=0; nEvent<nMC; ++nEvent){
    fitter_tree->GetEntry(nEvent);

    hMCnVtx->Fill(MCnVtx);

    hMCnVtx_aux->Fill(MCnVtx,normMC);

  }



  hWeightnVtx = hDatanVtx;
  hWeightnVtx->Divide(hMCnVtx);

  hWeightnVtx_aux = hDatanVtx_aux;
  hWeightnVtx_aux->Divide(hMCnVtx_aux);

  hDatanVtx->SetMarkerStyle(20);
  hDatanVtx2->SetMarkerStyle(20);
  hMCnVtx->SetLineColor(2);

  hDatanVtx_aux->SetMarkerStyle(20);
  hDatanVtx2_aux->SetMarkerStyle(20);
  hMCnVtx_aux->SetLineColor(2);

  hDatanVtx->GetXaxis()->SetTitle("# PV");
  hDatanVtx2->GetXaxis()->SetTitle("# PV");
  hMCnVtx->GetXaxis()->SetTitle("# PV");

  hDatanVtx->GetYaxis()->SetTitle("weight");
  hDatanVtx2->GetYaxis()->SetTitle("entries");
  hMCnVtx->GetYaxis()->SetTitle("entries");

  TCanvas *c = new TCanvas();
  c->Divide(2,2);
  c->cd(1);
  hDatanVtx2->Draw("pe");
  c->cd(2);
  hMCnVtx->Draw("");
  c->cd(3);
  hDatanVtx2->Draw("pe");
  hMCnVtx->SetLineColor(2);
  hMCnVtx->Draw("same");

  c->cd(4);
  hWeightnVtx->Draw("pe");


  TCanvas *ca = new TCanvas();
  ca->Divide(2,2);
  ca->cd(1);
  hDatanVtx2_aux->Draw("pe");
  ca->cd(2);
  hMCnVtx_aux->Draw("");
  ca->cd(3);
  hDatanVtx2_aux->Draw("pe");
  hMCnVtx_aux->SetLineColor(2);
  hMCnVtx_aux->Draw("same");

  ca->cd(4);
  hWeightnVtx_aux->Draw("pe");


  //inputMC->cd("tpTreeMuMu");
  inputMC->cd("tpTreeElElBDT");
  TBranch *newBranch = fitter_tree->Branch("weight", &weight, "weight/F");

  for(size_t nEvent=0; nEvent<nMC; ++nEvent){

    fitter_tree->GetEntry(nEvent);

    weight = hWeightnVtx_aux->GetBinContent(MCnVtx+1);
    hWeightedMC->Fill(MCnVtx,weight*normMC);
    newBranch->Fill();

  }

  TCanvas *cc = new TCanvas();
  cc->cd();
  hWeightedMC->SetLineColor(3);
  hWeightedMC->SetLineWidth(2);
  hWeightedMC->GetXaxis()->SetTitle("# PV");
  hWeightedMC->GetYaxis()->SetTitle("Weight");
  hWeightedMC->Draw("");
  hDatanVtx2_aux->Draw("pe same");

  inputMC->Write("", TObject::kOverwrite);

}

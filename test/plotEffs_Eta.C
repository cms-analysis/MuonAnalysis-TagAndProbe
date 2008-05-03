

{

   gStyle->SetFillColor(10);
   gStyle->SetTitleFillColor(10);
   gStyle->SetStatColor(10);
   gStyle->SetErrorX(0);
   gStyle->SetCanvasColor(10);

   TFile *_file0 = TFile::Open("muon_eff_test.root");
   _file0.cd();
   
   TH1F *heff_Eta = (TH1F*)_file0->Get("heff_Eta");
   TH1F *heff_sbs_Eta = (TH1F*)_file0->Get("heff_sbs_Eta");
   TH1F *truth_eff_Eta = (TH1F*)_file0->Get("truth_eff_Eta");

   heff_Eta->SetMinimum(0.9);
   heff_Eta->SetMaximum(1.05);
   heff_Eta->SetMarkerStyle(kFullCircle);
   heff_Eta->SetMarkerSize(0.5);
   heff_Eta->GetXaxis()->SetTitle("p_{T} (GeV/c);");
   heff_Eta->GetYaxis()->SetTitle("Efficiency");
   heff_Eta->SetTitle("Efficiency vs p_{T}");
   heff_Eta->SetStats(kFALSE);
   heff_Eta->Draw();

   heff_sbs_Eta->SetMarkerColor(kBlue);
   heff_sbs_Eta->SetMarkerStyle(kFullCircle);
   heff_sbs_Eta->SetMarkerSize(0.5);
   heff_sbs_Eta->SetLineColor(kBlue);
   heff_sbs_Eta->Draw("same");

   truth_eff_Eta->SetMarkerColor(kRed);
   truth_eff_Eta->SetMarkerStyle(kFullCircle);
   truth_eff_Eta->SetMarkerSize(0.5);
   truth_eff_Eta->SetLineColor(kRed);
   truth_eff_Eta->Draw("same");


   TLegend t(0.65,0.78,0.95,0.98);
   t.AddEntry("heff_Eta","Fit Result");
   t.AddEntry("heff_sbs_Eta","SB Sub Result");
   t.AddEntry("truth_eff_Eta","MC Truth");
   t.Draw();
}

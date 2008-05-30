

{

   gStyle->SetFillColor(10);
   gStyle->SetTitleFillColor(10);
   gStyle->SetStatColor(10);
   gStyle->SetErrorX(0);
   gStyle->SetCanvasColor(10);

   TFile *_file0 = TFile::Open("muon_eff_test.root");
   _file0.cd();
   
   TH1F *fit_eff_Eta = (TH1F*)_file0->Get("fit_eff_Eta");
   TH1F *sbs_eff_Eta = (TH1F*)_file0->Get("sbs_eff_Eta");
   TH1F *truth_eff_Eta = (TH1F*)_file0->Get("truth_eff_Eta");

   fit_eff_Eta->SetMinimum(0.9);
   fit_eff_Eta->SetMaximum(1.05);
   fit_eff_Eta->SetMarkerStyle(kFullCircle);
   fit_eff_Eta->SetMarkerSize(0.5);
   fit_eff_Eta->GetXaxis()->SetTitle("#eta");
   fit_eff_Eta->GetYaxis()->SetTitle("Efficiency");
   fit_eff_Eta->SetTitle("Efficiency vs Pseudorapidity");
   fit_eff_Eta->SetStats(kFALSE);
   fit_eff_Eta->Draw();

   sbs_eff_Eta->SetMarkerColor(kBlue);
   sbs_eff_Eta->SetMarkerStyle(kFullCircle);
   sbs_eff_Eta->SetMarkerSize(0.5);
   sbs_eff_Eta->SetLineColor(kBlue);
   sbs_eff_Eta->Draw("same");

   truth_eff_Eta->SetMarkerColor(kRed);
   truth_eff_Eta->SetMarkerStyle(kFullCircle);
   truth_eff_Eta->SetMarkerSize(0.5);
   truth_eff_Eta->SetLineColor(kRed);
   truth_eff_Eta->Draw("same");


   TLegend t(0.65,0.78,0.95,0.98);
   t.AddEntry("fit_eff_Eta","Fit Result");
   t.AddEntry("sbs_eff_Eta","SB Sub Result");
   t.AddEntry("truth_eff_Eta","MC Truth");
   t.Draw();
}

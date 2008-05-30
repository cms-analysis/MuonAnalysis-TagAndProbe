

{

   gStyle->SetFillColor(10);
   gStyle->SetTitleFillColor(10);
   gStyle->SetStatColor(10);
   gStyle->SetErrorX(0);
   gStyle->SetCanvasColor(10);

   TFile *_file0 = TFile::Open("muon_eff_test.root");
   _file0.cd();
   
   TH1F *fit_eff_Pt = (TH1F*)_file0->Get("fit_eff_Pt");
   TH1F *sbs_eff_Pt = (TH1F*)_file0->Get("sbs_eff_Pt");
   TH1F *truth_eff_Pt = (TH1F*)_file0->Get("truth_eff_Pt");

   fit_eff_Pt->SetMinimum(0.9);
   fit_eff_Pt->SetMaximum(1.05);
   fit_eff_Pt->SetMarkerStyle(kFullCircle);
   fit_eff_Pt->SetMarkerSize(0.5);
   fit_eff_Pt->GetXaxis()->SetTitle("p_{T} (GeV/c)");
   fit_eff_Pt->GetYaxis()->SetTitle("Efficiency");
   fit_eff_Pt->SetTitle("Efficiency vs p_{T}");
   fit_eff_Pt->SetStats(kFALSE);
   fit_eff_Pt->Draw();

   sbs_eff_Pt->SetMarkerColor(kBlue);
   sbs_eff_Pt->SetMarkerStyle(kFullCircle);
   sbs_eff_Pt->SetMarkerSize(0.5);
   sbs_eff_Pt->SetLineColor(kBlue);
   sbs_eff_Pt->Draw("same");

   truth_eff_Pt->SetMarkerColor(kRed);
   truth_eff_Pt->SetMarkerStyle(kFullCircle);
   truth_eff_Pt->SetMarkerSize(0.5);
   truth_eff_Pt->SetLineColor(kRed);
   truth_eff_Pt->Draw("same");

   TLegend t(0.65,0.78,0.95,0.98);
   t.AddEntry("fit_eff_Pt","Fit Result");
   t.AddEntry("sbs_eff_Pt","SB Sub Result");
   t.AddEntry("truth_eff_Pt","MC Truth");
   t.Draw();
}

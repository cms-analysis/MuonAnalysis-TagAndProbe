

{

   gStyle->SetFillColor(10);
   gStyle->SetTitleFillColor(10);
   gStyle->SetStatColor(10);
   gStyle->SetErrorX(0);
   gStyle->SetCanvasColor(10);

   TFile *_file0 = TFile::Open("muon_eff_test.root");
   _file0.cd();
   
   TH1F *heff_Pt = (TH1F*)_file0->Get("heff_Pt");
   TH1F *heff_sbs_Pt = (TH1F*)_file0->Get("heff_sbs_Pt");
   TH1F *truth_eff_Pt = (TH1F*)_file0->Get("truth_eff_Pt");

   heff_Pt->SetMinimum(0.9);
   heff_Pt->SetMaximum(1.05);
   heff_Pt->SetMarkerStyle(kFullCircle);
   heff_Pt->SetMarkerSize(0.5);
   heff_Pt->GetXaxis()->SetTitle("p_{T} (GeV/c);");
   heff_Pt->GetYaxis()->SetTitle("Efficiency");
   heff_Pt->SetTitle("Efficiency vs p_{T}");
   heff_Pt->SetStats(kFALSE);
   heff_Pt->Draw();

   heff_sbs_Pt->SetMarkerColor(kBlue);
   heff_sbs_Pt->SetMarkerStyle(kFullCircle);
   heff_sbs_Pt->SetMarkerSize(0.5);
   heff_sbs_Pt->SetLineColor(kBlue);
   heff_sbs_Pt->Draw("same");

   truth_eff_Pt->SetMarkerColor(kRed);
   truth_eff_Pt->SetMarkerStyle(kFullCircle);
   truth_eff_Pt->SetMarkerSize(0.5);
   truth_eff_Pt->SetLineColor(kRed);
   truth_eff_Pt->Draw("same");

   TLegend t(0.65,0.78,0.95,0.98);
   t.AddEntry("heff_Pt","Fit Result");
   t.AddEntry("heff_sbs_Pt","SB Sub Result");
   t.AddEntry("truth_eff_Pt","MC Truth");
   t.Draw();
}

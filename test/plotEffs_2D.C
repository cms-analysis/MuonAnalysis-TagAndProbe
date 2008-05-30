

{

   gStyle->SetFillColor(10);
   gStyle->SetTitleFillColor(10);
   gStyle->SetStatColor(10);
   gStyle->SetErrorX(0);
   gStyle->SetCanvasColor(10);

   TFile *_file0 = TFile::Open("muon_eff_test.root");
   _file0.cd();
   
   TH2F *fit_eff_Pt_Eta = (TH2F*)_file0->Get("fit_eff_Pt_Eta");
   TH2F *sbs_eff_Pt_Eta = (TH2F*)_file0->Get("sbs_eff_Pt_Eta");
   TH2F *truth_eff_Pt_Eta = (TH2F*)_file0->Get("truth_eff_Pt_Eta");

   // For each eta bin plot vs pt
   TCanvas *c[100];

   for( int ybin=1; ybin<=fit_eff_Pt_Eta->GetNbinsY(); ++ybin )
   {
      char cname[100];
      sprintf( cname, "cy%d", ybin );
      c[ybin-1] = new TCanvas(cname,"Plots",700,500);
      c[ybin-1]->SetFillColor(10);

      char hname[100];
      sprintf( hname, "hy_%d", ybin );
      TH1D *h = fit_eff_Pt_Eta->ProjectionX(hname,ybin,ybin);
      char hname_sbs[100];
      sprintf( hname_sbs, "hy_sbs_%d", ybin );
      TH1D *h_sbs = sbs_eff_Pt_Eta->ProjectionX(hname_sbs,ybin,ybin);
      char hname_truth[100];
      sprintf( hname_truth, "hy_truth_%d", ybin );
      TH1D *h_truth = truth_eff_Pt_Eta->ProjectionX(hname_truth,ybin,ybin);

      h->SetMinimum(0.9);
      h->SetMaximum(1.05);
      h->SetMarkerStyle(kFullCircle);
      h->SetMarkerSize(0.5);
      h->GetXaxis()->SetTitle("p_{T} (GeV/c)");
      h->GetYaxis()->SetTitle("Efficiency");
      char tname[200];
      sprintf( tname, "Efficiency vs p_{T}: #eta bin %d", ybin );
      h->SetTitle(tname);
      h->SetStats(kFALSE);
      h->Draw();

      h_sbs->SetMarkerColor(kBlue);
      h_sbs->SetMarkerStyle(kFullCircle);
      h_sbs->SetMarkerSize(0.5);
      h_sbs->SetLineColor(kBlue);
      h_sbs->Draw("same");

      h_truth->SetMarkerColor(kRed);
      h_truth->SetMarkerStyle(kFullCircle);
      h_truth->SetMarkerSize(0.5);
      h_truth->SetLineColor(kRed);
      h_truth->Draw("same");
      
      TLegend *t = new TLegend(0.65,0.78,0.95,0.98);
      t->AddEntry(hname,"Fit Result");
      t->AddEntry(hname_sbs,"SB Sub Result");
      t->AddEntry(hname_truth,"MC Truth");
      t->Draw();
   }

   TCanvas *cx[100];
   for( int xbin=1; xbin<=fit_eff_Pt_Eta->GetNbinsX(); ++xbin )
   {
      char cname[100];
      sprintf( cname, "cx%d", xbin );
      cx[xbin-1] = new TCanvas(cname,"Plots",700,500);
      cx[xbin-1]->SetFillColor(10);

      char hname[100];
      sprintf( hname, "hx_%d", xbin );
      TH1D *h = fit_eff_Pt_Eta->ProjectionY(hname,xbin,xbin);
      char hname_sbs[100];
      sprintf( hname_sbs, "hx_sbs_%d", xbin );
      TH1D *h_sbs = sbs_eff_Pt_Eta->ProjectionY(hname_sbs,xbin,xbin);
      char hname_truth[100];
      sprintf( hname_truth, "hx_truth_%d", xbin );
      TH1D *h_truth = truth_eff_Pt_Eta->ProjectionY(hname_truth,xbin,xbin);

      h->SetMinimum(0.9);
      h->SetMaximum(1.05);
      h->SetMarkerStyle(kFullCircle);
      h->SetMarkerSize(0.5);
      h->GetXaxis()->SetTitle("#eta");
      h->GetYaxis()->SetTitle("Efficiency");
      char tname[200];
      sprintf( tname, "Efficiency vs #eta: p_{T} bin %d", xbin );
      h->SetTitle(tname);
      h->SetStats(kFALSE);
      h->Draw();

      h_sbs->SetMarkerColor(kBlue);
      h_sbs->SetMarkerStyle(kFullCircle);
      h_sbs->SetMarkerSize(0.5);
      h_sbs->SetLineColor(kBlue);
      h_sbs->Draw("same");

      h_truth->SetMarkerColor(kRed);
      h_truth->SetMarkerStyle(kFullCircle);
      h_truth->SetMarkerSize(0.5);
      h_truth->SetLineColor(kRed);
      h_truth->Draw("same");
      
      TLegend *t = new TLegend(0.65,0.78,0.95,0.98);
      t->AddEntry(hname,"Fit Result");
      t->AddEntry(hname_sbs,"SB Sub Result");
      t->AddEntry(hname_truth,"MC Truth");
      t->Draw();
   }
}

TString prefix = "plots/tracking/signal_0.5pb/";

TDirectory *mc_pt_eta = gFile->GetDirectory("histoTracking/pt_eta_mcTrue/");
TDirectory *fit_pt_eta = gFile->GetDirectory("histoTracking/pt_eta/");
TDirectory *mc_eta_phi = gFile->GetDirectory("histoTracking/eta_phi_mcTrue/");
TDirectory *fit_eta_phi = gFile->GetDirectory("histoTracking/eta_phi/");

TCanvas *c1 = new TCanvas("c1","c1");
void plotTracking() {
    doCanvas(fit_pt_eta,  3, 2, "eta_pt_%d_%d",  "eta_bin%d__pt_bin%d__hasValidHits_pass__tag_HLTMu3_pass__gaussPlusCubic");
    doCanvas(fit_eta_phi, 3, 6, "eta_phi_%d_%d", "eta_bin%d__phi_bin%d__hasValidHits_pass__tag_HLTMu3_pass__gaussPlusCubic");

    gROOT->ProcessLine(".x ~/cpp/tdrstyle.cc");
    c1 = new TCanvas("c1","c1");


    stack(mc_pt_eta, fit_pt_eta, "pt_barrel", "pt_plot__eta_bin1__hasValidHits_pass__mcTrue_true__tag_HLTMu3_pass");
    stack(mc_pt_eta, fit_pt_eta, "pt_ec_neg", "pt_plot__eta_bin0__hasValidHits_pass__mcTrue_true__tag_HLTMu3_pass");
    stack(mc_pt_eta, fit_pt_eta, "pt_ec_pos", "pt_plot__eta_bin2__hasValidHits_pass__mcTrue_true__tag_HLTMu3_pass");
    stack(mc_pt_eta, fit_pt_eta, "eta_lowpt", "eta_plot__pt_bin0__hasValidHits_pass__mcTrue_true__tag_HLTMu3_pass");
    stack(mc_pt_eta, fit_pt_eta, "eta_higpt", "eta_plot__pt_bin0__hasValidHits_pass__mcTrue_true__tag_HLTMu3_pass");
    stack(mc_eta_phi, fit_eta_phi, "phi_barrel", "phi_plot__eta_bin1__hasValidHits_pass__mcTrue_true__tag_HLTMu3_pass");
    stack(mc_eta_phi, fit_eta_phi, "phi_ec_neg", "phi_plot__eta_bin0__hasValidHits_pass__mcTrue_true__tag_HLTMu3_pass");
    stack(mc_eta_phi, fit_eta_phi, "phi_ec_pos", "phi_plot__eta_bin2__hasValidHits_pass__mcTrue_true__tag_HLTMu3_pass");
}

void stack(TDirectory *mc, TDirectory *fit, TString alias, TString mcname) {
    c1->cd(); c1->Clear();
    RooPlot *pmc = mc->Get("cnt_eff_plots/"+mcname);
    RooHist *hmc = (RooHist *) pmc->findObject("hxy_cnt_eff");
    hmc->SetLineColor(kRed);
    hmc->SetMarkerColor(kRed);
    hmc->SetMarkerStyle(24);
    pmc->Draw();
    pmc->GetYaxis()->SetRangeUser(0.9, 1.02);

    TString fitname = mcname; fitname.ReplaceAll("__mcTrue_true","");
    RooPlot *pfit = (RooPlot *) fit->Get("fit_eff_plots/"+fitname);
    RooHist *hfit = (RooHist *) pfit->findObject("hxy_fit_eff");
    hfit->SetLineColor(kBlue);
    hfit->SetMarkerColor(kBlue);
    pfit->Draw("SAME");

    c1->Print(prefix+alias+".png");
}

void doCanvas(TDirectory *dir, int binsx, int binsy, const char * easyname, const char * truename) {
    char buff[1023], baff[1023];
    for (int i = 0; i < binsx; ++i) {
        for (int j = 0; j < binsy; ++j) {
            sprintf(buff,easyname,i,j);
            sprintf(baff,truename,i,j);
            TCanvas *fitc = (TCanvas *) dir->Get(TString(baff)+"/fit_canvas");
            fitc->Draw(); 
            fitc->Print(prefix+TString("canvases/")+buff+"_fit.png");
            TCanvas *distc = (TCanvas *) dir->Get(TString(baff)+"/distributions_canvas");
            distc->Draw(); 
            distc->Print(prefix+TString("canvases/")+buff+"_dist.png");
        }
    }
}

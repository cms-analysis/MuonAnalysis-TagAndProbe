TString prefix = "plots/tracking/all_0.1pb/";


TCanvas *c1 = new TCanvas("c1","c1");
void plotTracking(int mc=1) {
    if (mc) {
        plotTrackingMC();
    } else {
        plotTrackingData();
    }
}

void plotTrackingData() {
    prefix = "plots/tracking/data_0.001pb/";
    TDirectory *fit_pt_eta = gFile->GetDirectory("histoTracking/eff");

    doCanvas(fit_pt_eta,  2, 1, "eff_pt_%d",  "pt_bin%d__hasValidHits_pass__gaussPlusCubic");

    gROOT->ProcessLine(".x ~/cpp/tdrstyle.cc");
    c1 = new TCanvas("c1","c1");

    single(fit_pt_eta, "fit", "pt_plot__hasValidHits_pass");
}

void plotTrackingMC() {
    TDirectory *mc_pt_eta = gFile->GetDirectory("histoTracking/pt_eta_mcTrue/");
    TDirectory *fit_pt_eta = gFile->GetDirectory("histoTracking/pt_eta/");
    TDirectory *mc_eta_phi = gFile->GetDirectory("histoTracking/eta_phi_mcTrue/");
    TDirectory *fit_eta_phi = gFile->GetDirectory("histoTracking/eta_phi/");

    doCanvas(fit_pt_eta,  3, 2, "eta_pt_%d_%d",  "eta_bin%d__pt_bin%d__hasValidHits_pass__tag_HLTMu3_pass__gaussPlusCubic");
    doCanvas(fit_eta_phi, 3, 6, "eta_phi_%d_%d", "eta_bin%d__phi_bin%d__hasValidHits_pass__tag_HLTMu3_pass__gaussPlusCubic");

    gROOT->ProcessLine(".x ~/cpp/tdrstyle.cc");
    c1 = new TCanvas("c1","c1");

    stack(mc_pt_eta, fit_pt_eta, "pt_barrel", "pt_plot__eta_bin1__hasValidHits_pass__mcTrue_true__tag_HLTMu3_pass");
    stack(mc_pt_eta, fit_pt_eta, "pt_ec_neg", "pt_plot__eta_bin0__hasValidHits_pass__mcTrue_true__tag_HLTMu3_pass");
    stack(mc_pt_eta, fit_pt_eta, "pt_ec_pos", "pt_plot__eta_bin2__hasValidHits_pass__mcTrue_true__tag_HLTMu3_pass");
    stack(mc_pt_eta, fit_pt_eta, "eta_lowpt", "eta_plot__pt_bin0__hasValidHits_pass__mcTrue_true__tag_HLTMu3_pass");
    stack(mc_pt_eta, fit_pt_eta, "eta_higpt", "eta_plot__pt_bin1__hasValidHits_pass__mcTrue_true__tag_HLTMu3_pass");
    stack(mc_eta_phi, fit_eta_phi, "phi_barrel", "phi_plot__eta_bin1__hasValidHits_pass__mcTrue_true__tag_HLTMu3_pass");
    stack(mc_eta_phi, fit_eta_phi, "phi_ec_neg", "phi_plot__eta_bin0__hasValidHits_pass__mcTrue_true__tag_HLTMu3_pass");
    stack(mc_eta_phi, fit_eta_phi, "phi_ec_pos", "phi_plot__eta_bin2__hasValidHits_pass__mcTrue_true__tag_HLTMu3_pass");
}

void stack(TDirectory *mc, TDirectory *fit, TString alias, TString mcname) {
    c1->cd(); c1->Clear();
    RooPlot *pmc = mc->Get("cnt_eff_plots/"+mcname);
    RooHist *hmc = (RooHist *) pmc->findObject("hxy_cnt_eff");
    hmc->SetLineWidth(2);
    hmc->SetLineColor(kRed);
    hmc->SetMarkerColor(kRed);
    hmc->SetMarkerStyle(25);
    hmc->SetMarkerSize(2.0);
    pmc->Draw();
    pmc->GetYaxis()->SetRangeUser(0.9, 1.02);

    TString fitname = mcname; fitname.ReplaceAll("__mcTrue_true","");
    RooPlot *pfit = (RooPlot *) fit->Get("fit_eff_plots/"+fitname);
    RooHist *hfit = (RooHist *) pfit->findObject("hxy_fit_eff");
    hfit->SetLineWidth(2);
    hfit->SetLineColor(kBlack);
    hfit->SetMarkerColor(kBlack);
    hfit->SetMarkerStyle(20);
    hfit->SetMarkerSize(1.6);
    pfit->Draw("SAME");

    c1->Print(prefix+alias+".png");
}

void single( TDirectory *fit, TString alias, TString fitname) {
    RooPlot *pfit = (RooPlot *) fit->Get("fit_eff_plots/"+fitname);
    RooHist *hfit = (RooHist *) pfit->findObject("hxy_fit_eff");
    hfit->SetLineWidth(2);
    hfit->SetLineColor(kBlack);
    hfit->SetMarkerColor(kBlack);
    pfit->Draw("");
    c1->Print(prefix+alias+".png");
}

void doCanvas(TDirectory *dir, int binsx, int binsy, const char * easyname, const char * truename) {
    char buff[1023], baff[1023];
    for (int i = 0; i < binsx; ++i) {
        for (int j = 0; j < binsy; ++j) {
            if (binsx != 1 && binsy != 1) {
                sprintf(buff,easyname,i,j);
                sprintf(baff,truename,i,j);
            } else if (binsx != 1) {
                sprintf(buff,easyname,i);
                sprintf(baff,truename,i);
            } else if (binsy != 1) {
                sprintf(buff,easyname,j);
                sprintf(baff,truename,j);
            } else {
                sprintf(buff,easyname);
                sprintf(baff,truename);
            }
            TCanvas *fitc = (TCanvas *) dir->Get(TString(baff)+"/fit_canvas");
            fitc->Draw(); 
            fitc->Print(prefix+TString("canvases/")+buff+"_fit.png");
            TCanvas *distc = (TCanvas *) dir->Get(TString(baff)+"/distributions_canvas");
            distc->Draw(); 
            distc->Print(prefix+TString("canvases/")+buff+"_dist.png");
        }
    }
}

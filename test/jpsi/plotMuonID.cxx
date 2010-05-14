//TString prefix = "plots/muonid/signal_0.1pb/";
TString prefix = "FIXME";

TCanvas *c1 = new TCanvas("c1","c1");
void plotMuonID(TString scenario, int mc=1) {
    prefix = TString("plots/muonid/")+scenario+"/";
    gSystem->mkdir(prefix,true);

    gROOT->ProcessLine(".x ~/cpp/tdrstyle.cc");
    gStyle->SetOptStat(0);
    c1 = new TCanvas("c1","c1");

    if (mc) {
        plotMuonIDMC();
    } else {
        plotMuonIDData();
    }
}

void plotMuonIDData() {
    char *ids[3] = { "Glb", "Sta", "TMA" };
    for (size_t i = 0; i < 3; ++i) {
        TString idname(ids[i]);

        TDirectory *fit_pt     = gFile->GetDirectory("histoMuFromTk/POG_"+idname+"_pt/");
        TDirectory *fit_pt_eta = gFile->GetDirectory("histoMuFromTk/POG_"+idname+"_pt_eta/");
    
        single(fit_pt, idname+"_pt_all", "pt_plot__eta_bin0");

        single(fit_pt_eta, idname+"_pt_barrel", "pt_plot__eta_bin1");
        single(fit_pt_eta, idname+"_pt_ec_neg", "pt_plot__eta_bin0");
        single(fit_pt_eta, idname+"_pt_ec_pos", "pt_plot__eta_bin2");
        single(fit_pt_eta, idname+"_eta_pt2_4.5",  "eta_plot__pt_bin0");
        single(fit_pt_eta, idname+"_eta_pt4.5_15", "eta_plot__pt_bin1");

        doCanvas(fit_pt,     1, 2, idname+"_pt_%d",  "eta_bin0__pt_bin%d__gaussPlusExpo");
        doCanvas(fit_pt_eta, 3, 2, idname+"_eta_pt_%d_%d",  "eta_bin%d__pt_bin%d__gaussPlusExpo");
    }

}

void plotMuonIDMC() {
    char *ids[3] = { "Glb", "Sta", "TMA" };
    for (size_t i = 0; i < 3; ++i) {
        TString idname(ids[i]);

        TDirectory *mc_pt_eta  = gFile->GetDirectory("histoMuFromTk/POG_"+idname+"_pt_eta_mcTrue/");
        TDirectory *fit_pt_eta = gFile->GetDirectory("histoMuFromTk/POG_"+idname+"_pt_eta/");

        stack(mc_pt_eta, fit_pt_eta, idname+"_pt_barrel", "pt_plot__eta_bin1__mcTrue_true__tag_HLTMu3_pass");
        stack(mc_pt_eta, fit_pt_eta, idname+"_pt_ec_neg", "pt_plot__eta_bin0__mcTrue_true__tag_HLTMu3_pass");
        stack(mc_pt_eta, fit_pt_eta, idname+"_pt_ec_pos", "pt_plot__eta_bin2__mcTrue_true__tag_HLTMu3_pass");
        stack(mc_pt_eta, fit_pt_eta, idname+"_eta_pt2_3",   "eta_plot__pt_bin0__mcTrue_true__tag_HLTMu3_pass");
        stack(mc_pt_eta, fit_pt_eta, idname+"_eta_pt3_4.5", "eta_plot__pt_bin1__mcTrue_true__tag_HLTMu3_pass");
        stack(mc_pt_eta, fit_pt_eta, idname+"_eta_pt4.5_6", "eta_plot__pt_bin2__mcTrue_true__tag_HLTMu3_pass");
        stack(mc_pt_eta, fit_pt_eta, idname+"_eta_pt6_20",  "eta_plot__pt_bin3__mcTrue_true__tag_HLTMu3_pass");

        doCanvas(fit_pt_eta,  3, 4,  idname+"_eta_pt_%d_%d",  "eta_bin%d__pt_bin%d__tag_HLTMu3_pass__gaussPlusExpo");
    }
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
    pmc->GetYaxis()->SetRangeUser(0.0, 1.08);

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
    c1->cd(); c1->Clear();
    RooPlot *pfit = (RooPlot *) fit->Get("fit_eff_plots/"+fitname);
    RooHist *hfit = (RooHist *) pfit->findObject("hxy_fit_eff");
    hfit->SetLineWidth(2);
    hfit->SetLineColor(kBlack);
    hfit->SetMarkerColor(kBlack);
    hfit->SetMarkerStyle(20);
    hfit->SetMarkerSize(1.6);
    pfit->Draw();
    pfit->GetYaxis()->SetRangeUser(0.0, 1.08);

    c1->Print(prefix+alias+".png");
}

void prettyLine(TCanvas *canv, int pad, const char *cname, int color) {
    RooCurve *c = (RooCurve *) canv->GetPad(pad)->FindObject(cname);
    c->SetLineWidth(2);
    c->SetLineColor(color);
}
void prettyLines(TCanvas *c) {
   prettyLine(c, 1, "pdfPass_Norm[mass]",                      kRed  );
   prettyLine(c, 1, "pdfPass_Norm[mass]_Comp[backgroundPass]", kBlue );
   prettyLine(c, 2, "pdfFail_Norm[mass]",                      kRed  );
   prettyLine(c, 2, "pdfFail_Norm[mass]_Comp[backgroundFail]", kBlue );
   prettyLine(c, 3, "simPdf_Norm[mass]",                                     kRed  );
   prettyLine(c, 3, "simPdf_Norm[mass]_Comp[backgroundPass,backgroundFail]", kBlue );
}
void doCanvas(TDirectory *dir, int binsx, int binsy, const char * easyname, const char * truename) {
    gSystem->mkdir(prefix+"canvases/",true);
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
            if (fitc == 0) {
                std::cerr << "Didn't find " << TString(baff) << "/fit_canvas in " << dir->GetName() << std::endl;
                continue;
            }
            fitc->Draw(); 
            prettyLines(fitc);
            fitc->Print(prefix+TString("canvases/")+buff+"_fit.png");
            //TCanvas *distc = (TCanvas *) dir->Get(TString(baff)+"/distributions_canvas");
            //distc->Draw(); 
            //distc->Print(prefix+TString("canvases/")+buff+"_dist.png");
        }
    }
}

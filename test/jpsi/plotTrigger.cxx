//TString prefix = "plots/trigger/signal_0.1pb/";
TString prefix = "plots/trigger/all_0.1pb/";
bool includeCounting = true;

TCanvas *c1 = new TCanvas("c1","c1");
void plotTrigger(TString scenario, int mc=1) {
    prefix = TString("plots/trigger/")+scenario+"/";
    gSystem->mkdir(prefix,true);

    gROOT->ProcessLine(".x ~/cpp/tdrstyle.cc");
    gStyle->SetOptStat(0);
    c1 = new TCanvas("c1","c1");

    if (mc) {
        plotTriggerMC();
    } else {
        plotTriggerData();
    }
}

void plotTriggerData() {
    char *trigg[2] = { "HLTMu3", "L1DiMuOpen" };
    for (size_t i = 0; i < 2; ++i) {
        TString trigname(trigg[i]);

        TDirectory *fit_pt = gFile->GetDirectory("histoTrigger/"+trigname+"_pt/");
        TDirectory *fit_pt_eta = gFile->GetDirectory("histoTrigger/"+trigname+"_pt_eta/");
        TDirectory *fit_run = gFile->GetDirectory("histoTrigger/"+trigname+"_run/");
        TDirectory *fit_run_eta = gFile->GetDirectory("histoTrigger/"+trigname+"_run_eta/");

        single(fit_pt, trigname+"_pt_all", "pt_plot__eta_bin0__run_bin0__Glb_true");

        single(fit_pt_eta, trigname+"_pt_barrel", "pt_plot__eta_bin1__run_bin0__Glb_true");
        single(fit_pt_eta, trigname+"_pt_ec_neg", "pt_plot__eta_bin0__run_bin0__Glb_true");
        single(fit_pt_eta, trigname+"_pt_ec_pos", "pt_plot__eta_bin2__run_bin0__Glb_true");
        single(fit_pt_eta, trigname+"_eta_pt2_3",  "eta_plot__pt_bin0__run_bin0__Glb_true");
        single(fit_pt_eta, trigname+"_eta_pt3_12", "eta_plot__pt_bin1__run_bin0__Glb_true");

        single(fit_run, trigname+"_run_pt3", "run_plot__eta_bin0__pt_bin0__Glb_true");
        single(fit_run_eta, trigname+"_run_barrel_pt3", "run_plot__eta_bin1__pt_bin0__Glb_true");
        single(fit_run_eta, trigname+"_run_ec_neg_pt3", "run_plot__eta_bin0__pt_bin0__Glb_true");
        single(fit_run_eta, trigname+"_run_ec_pos_pt3", "run_plot__eta_bin2__pt_bin0__Glb_true");

        doCanvas(fit_pt,      1, 2, trigname+"_pt_%d",  "eta_bin0__pt_bin%d__run_bin0__Glb_true__gaussPlusExpo");
        doCanvas(fit_pt_eta,  3, 2, trigname+"_eta_pt_%d_%d",  "eta_bin%d__pt_bin%d__run_bin0__Glb_true__gaussPlusExpo");

    }

}

void plotTriggerMC() {
    char *trigg[2] = { "HLTMu3", "L1DiMuOpen" };
    for (size_t i = 0; i < 2; ++i) {
        TString trigname(trigg[i]);

        TDirectory *mc_pt_eta  = gFile->GetDirectory("histoTrigger/"+trigname+"_pt_eta_mcTrue/");
        TDirectory *fit_pt_eta = gFile->GetDirectory("histoTrigger/"+trigname+"_pt_eta/");

        stack(mc_pt_eta, fit_pt_eta, trigname+"_pt_barrel", "pt_plot__eta_bin1__Glb_true__mcTrue_true__tag_HLTMu3_pass");
        stack(mc_pt_eta, fit_pt_eta, trigname+"_pt_ec_neg", "pt_plot__eta_bin0__Glb_true__mcTrue_true__tag_HLTMu3_pass");
        stack(mc_pt_eta, fit_pt_eta, trigname+"_pt_ec_pos", "pt_plot__eta_bin2__Glb_true__mcTrue_true__tag_HLTMu3_pass");
        stack(mc_pt_eta, fit_pt_eta, trigname+"_eta_pt2_3",   "eta_plot__pt_bin0__Glb_true__mcTrue_true__tag_HLTMu3_pass");
        stack(mc_pt_eta, fit_pt_eta, trigname+"_eta_pt3_4.5", "eta_plot__pt_bin1__Glb_true__mcTrue_true__tag_HLTMu3_pass");
        stack(mc_pt_eta, fit_pt_eta, trigname+"_eta_pt4.5_6", "eta_plot__pt_bin2__Glb_true__mcTrue_true__tag_HLTMu3_pass");
        if (strstr(prefix.Data(), "signal_0.5pb")) {
            stack(mc_pt_eta, fit_pt_eta, trigname+"_eta_pt6_10",   "eta_plot__pt_bin3__Glb_true__mcTrue_true__tag_HLTMu3_pass");
            stack(mc_pt_eta, fit_pt_eta, trigname+"_eta_pt10_20",  "eta_plot__pt_bin4__Glb_true__mcTrue_true__tag_HLTMu3_pass");
            doCanvas(fit_pt_eta,  3, 5, trigname+"_eta_pt_%d_%d",  "eta_bin%d__pt_bin%d__Glb_true__tag_HLTMu3_pass__gaussPlusExpo");
        } else {
            stack(mc_pt_eta, fit_pt_eta, trigname+"_eta_pt6_12",  "eta_plot__pt_bin3__Glb_true__mcTrue_true__tag_HLTMu3_pass");
            doCanvas(fit_pt_eta,  3, 4, trigname+"_eta_pt_%d_%d",  "eta_bin%d__pt_bin%d__Glb_true__tag_HLTMu3_pass__gaussPlusExpo");
        }
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
    /*
    RooPlot *pcnt = (RooPlot *) fit->Get("cnt_eff_plots/"+fitname);
    RooHist *hcnt = (RooHist *) pcnt->findObject("hxy_cnt_eff");
    hcnt->SetLineWidth(2);
    hcnt->SetLineColor(kBlue);
    hcnt->SetMarkerColor(kBlue);
    hcnt->SetMarkerStyle(21);
    hcnt->SetMarkerSize(1.8);
    pcnt->Draw("SAME");
    */

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

    if (includeCounting) {
        RooPlot *pcnt = (RooPlot *) fit->Get("cnt_eff_plots/"+fitname);
        RooHist *hcnt = (RooHist *) pcnt->findObject("hxy_cnt_eff");
        hcnt->SetLineWidth(2);
        hcnt->SetLineColor(kBlue);
        hcnt->SetMarkerColor(kBlue);
        hcnt->SetMarkerStyle(21);
        hcnt->SetMarkerSize(1.8);
        pcnt->Draw();
        pcnt->GetYaxis()->SetRangeUser(0.0, 1.08);
    }

    RooPlot *pfit = (RooPlot *) fit->Get("fit_eff_plots/"+fitname);
    RooHist *hfit = (RooHist *) pfit->findObject("hxy_fit_eff");
    hfit->SetLineWidth(2);
    hfit->SetLineColor(kBlack);
    hfit->SetMarkerColor(kBlack);
    hfit->SetMarkerStyle(20);
    hfit->SetMarkerSize(1.6);
    if (includeCounting) {
        pfit->Draw("SAME");
    } else {
        pfit->Draw();
        pfit->GetYaxis()->SetRangeUser(0.0, 1.08);
    }

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

TString prefix = "plots/trigger/signal_0.1pb/";

TCanvas *c1 = new TCanvas("c1","c1");
void plotTrigger(int mc=1) {
    gSystem->mkdir(prefix,true);
    if (mc) {
        plotTriggerMC();
    } else {
        plotTriggerData();
    }
}

void plotTriggerData() {
    prefix = "plots/trigger/data_0.001pb/";
    gSystem->mkdir(prefix,true);
    char *trigg[2] = { "HLTMu3", "L1DiMuOpen" };

    for (size_t i = 0; i < 2; ++i) {
        TString trigname(trigg[i]);
        TDirectory *fit_pt_eta = gFile->GetDirectory("histoTrigger/"+trigname+"_pt_eta/");
        doCanvas(fit_pt_eta,  3, 2, trigname+"_eta_pt_%d_%d",  "eta_bin%d__pt_bin%d__run_bin0__Glb_true__gaussPlusExpo");

        TDirectory *fit_pt = gFile->GetDirectory("histoTrigger/"+trigname+"_pt/");
        doCanvas(fit_pt,  1, 2, trigname+"_pt_%d",  "eta_bin0__pt_bin%d__run_bin0__Glb_true__gaussPlusExpo");
    }

    gROOT->ProcessLine(".x ~/cpp/tdrstyle.cc");
    c1 = new TCanvas("c1","c1");

    for (size_t i = 0; i < 2; ++i) {
        TString trigname(trigg[i]);
        TDirectory *fit_pt_eta = gFile->GetDirectory("histoTrigger/"+trigname+"_pt_eta/");

        single(fit_pt_eta, trigname+"_pt_barrel", "pt_plot__eta_bin1__run_bin0__Glb_true");
        single(fit_pt_eta, trigname+"_pt_ec_neg", "pt_plot__eta_bin0__run_bin0__Glb_true");
        single(fit_pt_eta, trigname+"_pt_ec_pos", "pt_plot__eta_bin2__run_bin0__Glb_true");
        single(fit_pt_eta, trigname+"_eta_pt2_3",  "eta_plot__pt_bin0__run_bin0__Glb_true");
        single(fit_pt_eta, trigname+"_eta_pt3_15", "eta_plot__pt_bin1__run_bin0__Glb_true");

        TDirectory *fit_pt = gFile->GetDirectory("histoTrigger/"+trigname+"_pt/");
        single(fit_pt, trigname+"_pt_all", "pt_plot__eta_bin0__run_bin0__Glb_true");

    }

}

void plotTriggerMC() {
    //doCanvas(fit_pt_eta,  3, 2, "eta_pt_%d_%d",  "eta_bin%d__pt_bin%d__Glb_true__tag_HLTMu3_pass__gaussPlusCubic");

    gROOT->ProcessLine(".x ~/cpp/tdrstyle.cc");
    c1 = new TCanvas("c1","c1");

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
        stack(mc_pt_eta, fit_pt_eta, trigname+"_eta_pt6_20",  "eta_plot__pt_bin3__Glb_true__mcTrue_true__tag_HLTMu3_pass");
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
    RooPlot *pcnt = (RooPlot *) fit->Get("cnt_eff_plots/"+fitname);
    RooHist *hcnt = (RooHist *) pcnt->findObject("hxy_cnt_eff");
    hcnt->SetLineWidth(2);
    hcnt->SetLineColor(kBlue);
    hcnt->SetMarkerColor(kBlue);
    hcnt->SetMarkerStyle(21);
    hcnt->SetMarkerSize(1.8);
    pcnt->Draw("SAME");

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
    RooPlot *pcnt = (RooPlot *) fit->Get("cnt_eff_plots/"+fitname);
    RooHist *hcnt = (RooHist *) pcnt->findObject("hxy_cnt_eff");
    hcnt->SetLineWidth(2);
    hcnt->SetLineColor(kBlue);
    hcnt->SetMarkerColor(kBlue);
    hcnt->SetMarkerStyle(21);
    hcnt->SetMarkerSize(1.8);
    pcnt->Draw();
    pcnt->GetYaxis()->SetRangeUser(0.0, 1.08);

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
            fitc->Print(prefix+TString("canvases/")+buff+"_fit.png");
            //TCanvas *distc = (TCanvas *) dir->Get(TString(baff)+"/distributions_canvas");
            //distc->Draw(); 
            //distc->Print(prefix+TString("canvases/")+buff+"_dist.png");
        }
    }
}

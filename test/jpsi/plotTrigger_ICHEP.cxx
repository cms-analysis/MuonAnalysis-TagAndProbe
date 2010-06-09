/**
  USAGE: 
    DATA-only:    root.exe -b -l -q TnP_ICHEP_Trigger_data_all.root "plotTrigger.cxx(\"data_all\")"
    MC data-like: root.exe -b -l -q TnP_ICHEP_Trigger_datalike_mc.root  "plotTrigger.cxx(\"datalike_mc\")"
    DATA+MC:      root.exe -b -l -q TnP_ICHEP_Trigger_data_all.root TnP_ICHEP_Trigger_datalike_mc.root  "plotTrigger.cxx(\"data_vs_mc\")"

  REQUIRES:
   1) mkdir -p plots_ichep/muonid/
   2) provide a suitable "tdrStyle.cc" macro or similar
      (by default, it's taken from ~/cpp/tdrstyle.cc;
       if you need one you might want to grab ~gpetrucc/cpp/tdrstyle.cc)
*/
#include <TCanvas.h>
#include <TPad.h>
#include "plotUtil.cxx"
TString prefix = "plots_ichep_dev/trigger/";
TString basedir = "histoTrigger";

TFile *ref = 0;

TCanvas *c1 = 0;
void plotTrigger_ICHEP(TString scenario) {
    prefix = prefix+scenario+"/";
    gSystem->mkdir(prefix,true);

    gROOT->ProcessLine(".x /afs/cern.ch/user/g/gpetrucc/cpp/tdrstyle.cc");
    gStyle->SetOptStat(0);
    c1 = new TCanvas("c1","c1");

    if (gROOT->GetListOfFiles()->GetEntries() == 2) {
        ref = (TFile *) gROOT->GetListOfFiles()->At(1);
        ((TFile*) gROOT->GetListOfFiles()->At(0))->cd();
    }
    plotTriggerData();
}

void plotTriggerData() {
    char *ids[2] = { "Mu3", "L1DoubleMuOpen" };
    for (size_t i = 0; i < 2; ++i) {
        TString idname(ids[i]);
        plotTriggerData(idname);

    }

    for (size_t i = 0; i < 1; ++i) {
        TString idname = TString(ids[i])+"overL1";
        plotTriggerData(idname);
    }
    
}
void plotTriggerData(TString idname) {
    TDirectory *fit_pt_eta = gFile->GetDirectory(basedir+"/"+idname+"_pt_abseta/");
    if (ref != 0) {
        TDirectory *ref_pt_eta = ref->GetDirectory(basedir+"/"+idname+"_pt_abseta/");
        refstack(fit_pt_eta, ref_pt_eta, idname+"_pt_eta-0.0-0.8",  "pt_PLOT_abseta_bin0_");
        refstack(fit_pt_eta, ref_pt_eta, idname+"_pt_eta-0.8-1.2",  "pt_PLOT_abseta_bin1_");
        refstack(fit_pt_eta, ref_pt_eta, idname+"_pt_eta-1.2-1.6",  "pt_PLOT_abseta_bin2_");
        refstack(fit_pt_eta, ref_pt_eta, idname+"_pt_eta-1.6-2.4",  "pt_PLOT_abseta_bin3_");
        refstack(fit_pt_eta, ref_pt_eta, idname+"_eta_pt-3-5",      "abseta_PLOT_pt_bin1_");
        refstack(fit_pt_eta, ref_pt_eta, idname+"_eta_pt-5-12",     "abseta_PLOT_pt_bin2_");
    } else {
        TDirectory *mc_pt_eta  = gFile->GetDirectory(basedir+"/"+idname+"_pt_abseta_mcTrue/");
        if (mc_pt_eta) {
            mcstack(fit_pt_eta, mc_pt_eta, idname+"_pt_eta-0.0-0.8", "pt_PLOT_abseta_bin0_");
            mcstack(fit_pt_eta, mc_pt_eta, idname+"_pt_eta-0.8-1.2", "pt_PLOT_abseta_bin1_");
            mcstack(fit_pt_eta, mc_pt_eta, idname+"_pt_eta-1.2-1.6", "pt_PLOT_abseta_bin2_");
            mcstack(fit_pt_eta, mc_pt_eta, idname+"_pt_eta-1.6-2.4", "pt_PLOT_abseta_bin3_");
            mcstack(fit_pt_eta, mc_pt_eta, idname+"_eta_pt-3-5",     "abseta_PLOT_pt_bin1_");
            mcstack(fit_pt_eta, mc_pt_eta, idname+"_eta_pt-5-12",    "abseta_PLOT_pt_bin2_");
        } else {
            single(fit_pt_eta, idname+"_pt_eta-0.0-0.8", "pt_PLOT_abseta_bin0_");
            single(fit_pt_eta, idname+"_pt_eta-0.8-1.2", "pt_PLOT_abseta_bin1_");
            single(fit_pt_eta, idname+"_pt_eta-1.2-1.6", "pt_PLOT_abseta_bin2_");
            single(fit_pt_eta, idname+"_pt_eta-1.6-2.4", "pt_PLOT_abseta_bin3_");
            single(fit_pt_eta, idname+"_eta_pt-3-5",     "abseta_PLOT_pt_bin1_");
            single(fit_pt_eta, idname+"_eta_pt-5-12",    "abseta_PLOT_pt_bin2_");
        }
    }

    doCanvas(fit_pt_eta, 4, 3, idname+"_eta_%d_pt_%d",   "abseta_bin%d__pt_bin%d_");
}

void refstack(TDirectory *fit, TDirectory *ref, TString alias, TString fitname) {
    RooPlot *pref = getFromPrefix(ref->GetDirectory("fit_eff_plots"), fitname);
    if (pref == 0) {
        std::cerr << "NOT FOUND: " << "fit_eff_plots/"+fitname << " in " << ref->GetName() << std::endl;
        return;
    }
    RooHist *href = (RooHist *) pref->FindObject("hxy_fit_eff");
    href->SetLineWidth(2);
    href->SetLineColor(kRed);
    href->SetMarkerColor(kRed);
    href->SetMarkerStyle(25);
    href->SetMarkerSize(2.0);

    RooPlot *pfit = getFromPrefix(fit->GetDirectory("fit_eff_plots"), fitname);
    if (pfit == 0) {
        std::cerr << "NOT FOUND: " << "fit_eff_plots/"+fitname << " in " << fit->GetName() << std::endl;
        return;
    }
    RooHist *hfit = (RooHist *) pfit->FindObject("hxy_fit_eff");
    hfit->SetLineWidth(2);
    hfit->SetLineColor(kBlack);
    hfit->SetMarkerColor(kBlack);
    hfit->SetMarkerStyle(20);
    hfit->SetMarkerSize(1.6);

    pref->Draw();
    hfit->Draw("P SAME");
    gPad->Print(prefix+alias+".png");
}

void mcstack(TDirectory *fit, TDirectory *ref, TString alias, TString mcname) {
    RooPlot *pref = getFromPrefix(ref->GetDirectory("cnt_eff_plots"), mcname);
    RooHist *href = (RooHist *) pref->FindObject("hxy_cnt_eff");
    href->SetLineWidth(2);
    href->SetLineColor(209);
    href->SetMarkerColor(209);
    href->SetMarkerStyle(25);
    href->SetMarkerSize(2.0);

    RooPlot *pfit = getFromPrefix(fit->GetDirectory("fit_eff_plots"), fitname);
    RooHist *hfit = (RooHist *) pfit->FindObject("hxy_fit_eff");
    hfit->SetLineWidth(2);
    hfit->SetLineColor(206);
    hfit->SetMarkerColor(206);
    hfit->SetMarkerStyle(20);
    hfit->SetMarkerSize(1.6);

    pref->Draw();
    hfit->Draw("P SAME");
    gPad->Print(prefix+alias+".png");
}



void single( TDirectory *fit, TString alias, TString fitname) {
    RooPlot *pfit = getFromPrefix(fit->GetDirectory("fit_eff_plots"), fitname);
    if (pfit == 0) {
        std::cerr << "NOT FOUND: " << "fit_eff_plots/"+fitname << " in " << fit->GetName() << std::endl;
        return;
    }
    RooHist *hfit = (RooHist *) pfit->FindObject("hxy_fit_eff");
    hfit->SetLineWidth(2);
    hfit->SetLineColor(kBlue);
    hfit->SetMarkerColor(kBlue);
    hfit->SetMarkerStyle(20);
    hfit->SetMarkerSize(1.6);
    pfit->Print(prefix+alias+".png"); 
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
            TDirectory *subdir = (TDirectory *) getFromPrefix(dir, baff);
            if (subdir == 0) {
                std::cerr << "Didn't find '" << baff << "*' in " << dir->GetName() << std::endl;
                continue;
            }
            TCanvas *fitc = (TCanvas *) subdir->Get("fit_canvas");
            if (fitc == 0) {
                std::cerr << "Didn't find " << TString(baff) << "/fit_canvas in " << dir->GetName() << std::endl;
                continue;
            }
            fitc->Draw(); 
            prettyLines(fitc);
            fitc->Print(prefix+TString("canvases/")+buff+"_fit.png");
        }
    }
}

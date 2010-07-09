/**
  USAGE: 
    DATA-only:    root.exe -b -l -q TnP_ICHEP_Trigger_data_all.root 'plotTrigger_ICHEP.cxx("data_all")'
    MC data-like: root.exe -b -l -q TnP_ICHEP_Trigger_datalike_mc.root  'plotTrigger_ICHEP.cxx("datalike_mc")'
    DATA+MC:      root.exe -b -l -q TnP_ICHEP_Trigger_data_all.root TnP_ICHEP_Trigger_datalike_mc.root  'plotTrigger_ICHEP.cxx("data_vs_mc")'
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

    doRatioPlot = false;
    doDiffPlot = false;
    doPdf = true; 
    doSquare = true; yMax = 1.1;
    datalbl = "Data, 77 nb^{-1}";
    reflbl  = "Simulation";
    preliminary = "CMS Preliminary,   #sqrt{s} = 7 TeV";

    plotTriggerData();
}

bool extended = false;
void plotTriggerData() {
    TString mu[2] = { "Cal", "VBTFLike" };
    TString trig[2] = { "Mu3", "L1DoubleMuOpen" };
    for (size_t i = 0; i < 2; ++i) {
      for (size_t j = 0; j < 2; ++j) {
        if (i == 1 && j > 0 && !extended) continue;
        retitle = TString(j == 0 ? "HLT+L1" : "L1")+" efficiency";
        TString idname = mu[i]+"_To_"+trig[j];
        plotTriggerData(idname);
      }
      if (extended) {
          for (size_t j = 0; j < 1; ++j) {
            retitle = "HLT only efficiency";
            TString idname = mu[i]+"_To_"+trig[j]+"overL1";
            plotTriggerData(idname);
          }
      }
    }
}
void plotTriggerData(TString idname) {
    bool mu3 = (strstr(idname.Data(), "Mu3"));
    TDirectory *fit_pt  = gFile->GetDirectory(basedir+"/"+idname+"_pt/");
    TDirectory *fit_eta = gFile->GetDirectory(basedir+"/"+idname+"_abseta/");
    TDirectory *fit_c_r  = gFile->GetDirectory(basedir+"/"+idname+"_check_run/");
    TDirectory *fit_c_e  = gFile->GetDirectory(basedir+"/"+idname+"_check_eta/");
    if (ref != 0) {
        TDirectory *ref_pt  = ref->GetDirectory(basedir+"/"+idname+"_pt/");
        TDirectory *ref_eta = ref->GetDirectory(basedir+"/"+idname+"_abseta/");
        TDirectory *mc_pt   = ref->GetDirectory(basedir+"/"+idname+"_pt_mcTrue/");
        TDirectory *mc_eta  = ref->GetDirectory(basedir+"/"+idname+"_abseta_mcTrue/");

        extraSpam = "        |#eta| < 1.2"; refstack(fit_pt,  ref_pt,  idname+"_pt_barrel",   "pt_PLOT_abseta_bin0_");
        extraSpam = "  1.2 < |#eta| < 2.1"; refstack(fit_pt,  ref_pt,  idname+"_pt_endcaps",  "pt_PLOT_abseta_bin1_");
        if (extended && fit_eta) { extraSpam = "p_{T} = 3-5 GeV/c";  refstack(fit_eta, ref_eta, idname+"_eta_pt-3-5",  "abseta_PLOT_pt_bin1_"); }
        if (extended && fit_eta) { extraSpam = "p_{T} = 5-20 GeV/c"; refstack(fit_eta, ref_eta, idname+"_eta_pt-5-20", "abseta_PLOT_pt_bin2_"); }
        if (mc_pt) {
            extraSpam = "            |#eta| < 1.2"; refstack3(fit_pt, ref_pt, mc_pt, idname+"_pt_barrel_3",  "pt_PLOT_abseta_bin0_");
            extraSpam = "       1.2 < |#eta| < 2.1"; refstack3(fit_pt, ref_pt, mc_pt, idname+"_pt_endcaps_3", "pt_PLOT_abseta_bin1_");
        } 
    } else {
        TDirectory *mc_pt  = gFile->GetDirectory(basedir+"/"+idname+"_pt_mcTrue/");
        TDirectory *mc_eta = gFile->GetDirectory(basedir+"/"+idname+"_abseta_mcTrue/");
        if (mc_pt) {
            datalbl = "T&P fit"; reflbl = "Sim. truth";
            extraSpam = "        |#eta| < 1.2"; mcstack(fit_pt,  mc_pt,  idname+"_pt_barrel",   "pt_PLOT_abseta_bin0_");
            extraSpam = "  1.2 < |#eta| < 2.1"; mcstack(fit_pt,  mc_pt,  idname+"_pt_endcaps",  "pt_PLOT_abseta_bin1_");
            //if (!mu3) mcstack(fit_eta, mc_eta, idname+"_eta_pt-2-3",  "abseta_PLOT_pt_bin0_");
            if (extended && fit_eta) { extraSpam = "p_{T} = 3-5 GeV/c";  mcstack(fit_eta, mc_eta, idname+"_eta_pt-3-5",  "abseta_PLOT_pt_bin1_"); }
            if (extended && fit_eta) { extraSpam = "p_{T} = 5-20 GeV/c"; mcstack(fit_eta, mc_eta, idname+"_eta_pt-5-20", "abseta_PLOT_pt_bin2_"); }
        } else {
            single(fit_pt,  idname+"_pt_barrel",   "pt_PLOT_abseta_bin0_");
            single(fit_pt,  idname+"_pt_endcaps",  "pt_PLOT_abseta_bin1_");
            //if (!mu3) single(fit_eta, idname+"_eta_pt-2-3",  "abseta_PLOT_pt_bin0_");
            if (extended && fit_eta) single(fit_eta, idname+"_eta_pt-3-5",  "abseta_PLOT_pt_bin1_");
            if (extended && fit_eta) single(fit_eta, idname+"_eta_pt-5-20", "abseta_PLOT_pt_bin2_");
        }
    }
    
    preliminary = "Cross check 4 < pt < 7";
    if (extended && fit_c_e != 0) single(fit_c_e,  idname+"_check_eta",   "eta_PLOT_");
    preliminary = "Cross check 1.2 < |#eta| < 2.1, 4 < pt < 7";
    if (extended && fit_c_r != 0) single(fit_c_r,  idname+"_check_run",   "run_PLOT_");
    preliminary = "CMS Preliminary,   #sqrt{s} = 7 TeV";

    if (ref == 0) {
        if (extended && fit_eta) doCanvas(fit_eta, 5, 3, idname+"_ETA_eta_%d_pt_%d",  "abseta_bin%d__pt_bin%d_");
        doCanvas(fit_pt,  2, 5, idname+"_PT_eta_%d_pt_%d",   "abseta_bin%d__pt_bin%d_");
    }
}



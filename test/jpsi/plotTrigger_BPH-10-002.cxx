/**
  USAGE: 
    DATA-only:    root.exe -b -l -q TnP_ICHEP_Trigger_data_all.root 'plotTrigger_ICHEP.cxx("data_all")'
    MC data-like: root.exe -b -l -q TnP_ICHEP_Trigger_datalike_mc.root  'plotTrigger_ICHEP.cxx("datalike_mc")'
    DATA+MC:      root.exe -b -l -q TnP_ICHEP_Trigger_data_all.root TnP_ICHEP_Trigger_datalike_mc.root  'plotTrigger_ICHEP.cxx("data_vs_mc")'
*/
#include <TCanvas.h>
#include <TPad.h>
#include "plotUtil.cxx"
TString prefix = "plots_bph_10_002_dev/trigger/";
TString basedir = "histoTrigger";

TFile *ref = 0;

TCanvas *c1 = 0;
void plotTrigger_BPH_10_002(TString scenario) {
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
    //TString mu[4] = { "POG_Glb", "POG_GlbPT", "POG_TMA", "POG_TMLSAT" };
    TString mu[3] = { "Glb", "TM_Incl", "TM_Excl" };
    TString trig[2] = { "Mu3", "L1DoubleMuOpen" };
    for (size_t i = 0; i < 3; ++i) {
      for (size_t j = 0; j < 2; ++j) {
        TString idname = mu[i]+"_To_"+trig[j];
        plotTriggerData(idname);
      }
    }
}
void plotTriggerData(TString idname) {
    TDirectory *fit_pt_eta = gFile->GetDirectory(basedir+"/"+idname+"_pt_abseta/");
    if (ref != 0) {
        TDirectory *ref_pt_eta = ref->GetDirectory(basedir+"/"+idname+"_pt_abseta/");
        refstack(fit_pt_eta, ref_pt_eta, idname+"_pt_barrel",  "pt_PLOT_abseta_bin0_");
        refstack(fit_pt_eta, ref_pt_eta, idname+"_pt_endcaps", "pt_PLOT_abseta_bin1_");
    } else {
        TDirectory *mc_pt_eta  = gFile->GetDirectory(basedir+"/"+idname+"_pt_abseta_mcTrue/");
        if (mc_pt_eta) {
            mcstack(fit_pt_eta, mc_pt_eta, idname+"_pt_barrel",  "pt_PLOT_abseta_bin0_");
            mcstack(fit_pt_eta, mc_pt_eta, idname+"_pt_endcaps", "pt_PLOT_abseta_bin1_");
        } else {
            single(fit_pt_eta, idname+"_pt_barrel",  "pt_PLOT_abseta_bin0_");
            single(fit_pt_eta, idname+"_pt_endcaps", "pt_PLOT_abseta_bin1_");
        }
    }

    doCanvas(fit_pt_eta, 1, 5, idname+"_barrel_pt_%d",   "abseta_bin0__pt_bin%d_");
    doCanvas(fit_pt_eta, 1, 5, idname+"_endcaps_pt_%d",  "abseta_bin1__pt_bin%d_");
}



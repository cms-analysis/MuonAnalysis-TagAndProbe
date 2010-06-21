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
    plotTriggerData();
}

void plotTriggerData() {
    //TString mu[4] = { "POG_Glb", "POG_GlbPT", "POG_TMA", "POG_TMLSAT" };
    TString mu[3] = { "POG_Glb", "Cal", "VBTFLike" };
    //TString mu[1] = { "POG_Glb" };
    TString trig[2] = { "Mu3", "L1DoubleMuOpen" };
    for (size_t i = 0; i < 3; ++i) {
      for (size_t j = 0; j < 2; ++j) {
        TString idname = mu[i]+"_To_"+trig[j];
        plotTriggerData(idname);
      }

      for (size_t j = 0; j < 1; ++j) {
        TString idname = mu[i]+"_To_"+trig[j]+"overL1";
        plotTriggerData(idname);
      }
    }
}
void plotTriggerData(TString idname) {
    bool mu3 = (strstr(idname.Data(), "Mu3"));
    TDirectory *fit_pt  = gFile->GetDirectory(basedir+"/"+idname+"_pt/");
    TDirectory *fit_eta = gFile->GetDirectory(basedir+"/"+idname+"_abseta/");
    if (ref != 0) {
        TDirectory *ref_pt  = ref->GetDirectory(basedir+"/"+idname+"_pt/");
        TDirectory *ref_eta = ref->GetDirectory(basedir+"/"+idname+"_abseta/");

        refstack(fit_pt,  ref_pt,  idname+"_pt_barrel",   "pt_PLOT_abseta_bin0_");
        refstack(fit_pt,  ref_pt,  idname+"_pt_endcaps",  "pt_PLOT_abseta_bin1_");
        //if (!mu3) refstack(fit_eta, ref_eta, idname+"_eta_pt-2-3",  "abseta_PLOT_pt_bin0_");
        refstack(fit_eta, ref_eta, idname+"_eta_pt-3-5",  "abseta_PLOT_pt_bin1_");
        refstack(fit_eta, ref_eta, idname+"_eta_pt-5-20", "abseta_PLOT_pt_bin2_");
    } else {
        TDirectory *mc_pt  = gFile->GetDirectory(basedir+"/"+idname+"_pt_mcTrue/");
        TDirectory *mc_eta = gFile->GetDirectory(basedir+"/"+idname+"_abseta_mcTrue/");
        if (mc_pt) {
            mcstack(fit_pt,  mc_pt,  idname+"_pt_barrel",   "pt_PLOT_abseta_bin0_");
            mcstack(fit_pt,  mc_pt,  idname+"_pt_endcaps",  "pt_PLOT_abseta_bin1_");
            //if (!mu3) mcstack(fit_eta, mc_eta, idname+"_eta_pt-2-3",  "abseta_PLOT_pt_bin0_");
            mcstack(fit_eta, mc_eta, idname+"_eta_pt-3-5",  "abseta_PLOT_pt_bin1_");
            mcstack(fit_eta, mc_eta, idname+"_eta_pt-5-20", "abseta_PLOT_pt_bin2_");
        } else {
            single(fit_pt,  idname+"_pt_barrel",   "pt_PLOT_abseta_bin0_");
            single(fit_pt,  idname+"_pt_endcaps",  "pt_PLOT_abseta_bin1_");
            //if (!mu3) single(fit_eta, idname+"_eta_pt-2-3",  "abseta_PLOT_pt_bin0_");
            single(fit_eta, idname+"_eta_pt-3-5",  "abseta_PLOT_pt_bin1_");
            single(fit_eta, idname+"_eta_pt-5-20", "abseta_PLOT_pt_bin2_");
        }
    }

    if (ref == 0) {
        doCanvas(fit_eta, 5, 3, idname+"_ETA_eta_%d_pt_%d",  "abseta_bin%d__pt_bin%d_");
        doCanvas(fit_pt,  2, 4, idname+"_PT_eta_%d_pt_%d",   "abseta_bin%d__pt_bin%d_");
    }
}



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

    datalbl = "Data, 84 nb^{-1}";
    preliminary = "";
    doRatioPlot = false; doDiffPlot = false;
    doPdf = false;
    doLogX = true;
    plotTriggerData();
}

void plotTriggerData() {
    TString mu[3] = { "Glb", "TM_Incl", "TM_Excl" };
    TString trig[2] = { "L1DoubleMuOpen", "Mu3"};
    for (size_t i = 0; i < 3; ++i) {
      for (size_t j = 0; j < 1; ++j) {
        TString idname = mu[i]+"_To_"+trig[j];
        plotTriggerData(idname);
      }
    }
}
void plotTriggerData(TString idname) {
    TDirectory *fit_pt = gFile->GetDirectory(basedir+"/"+idname+"_pt_abseta/");
    bool abseta = true;
    if (ref != 0) {
        TDirectory *ref_pt = ref->GetDirectory(basedir+"/"+idname+"_pt_abseta/");
        if (abseta) {
            refstack(fit_pt, ref_pt, idname+"_pt_barrel",  "pt_PLOT_abseta_bin0_");
            refstack(fit_pt, ref_pt, idname+"_pt_endcaps", "pt_PLOT_abseta_bin1_");
        } else {
            refstack(fit_pt, ref_pt, idname+"_pt_endcaps_neg", "pt_PLOT_eta_bin0_");
            refstack(fit_pt, ref_pt, idname+"_pt_barrel_neg",  "pt_PLOT_eta_bin1_");
            refstack(fit_pt, ref_pt, idname+"_pt_barrel_pos",  "pt_PLOT_eta_bin2_");
            refstack(fit_pt, ref_pt, idname+"_pt_endcaps_pos", "pt_PLOT_eta_bin3_");
        }
    } else {
        TDirectory *mc_pt = gFile->GetDirectory(basedir+"/"+idname+"_pt_abseta_mcTrue/");
        if (mc_pt) {
            mcstack(fit_pt, mc_pt, idname+"_pt_barrel",  "pt_PLOT_abseta_bin0_");
            mcstack(fit_pt, mc_pt, idname+"_pt_endcaps", "pt_PLOT_abseta_bin1_");
        } else {
            single(fit_pt, idname+"_pt_barrel",  "pt_PLOT_abseta_bin0_");
            single(fit_pt, idname+"_pt_endcaps", "pt_PLOT_abseta_bin1_");
        }
    }

    TDirectory *fit_vtx = gFile->GetDirectory(basedir+"/"+idname+"_vtx/");
    if ((ref == 0) && (fit_vtx != 0)) {
        doLogX = false;
        char buff[1255], baff[1255];
        for (int b = 0; b < 4; ++b) { for (int be = 0; be <= 1; ++be) {
            sprintf(buff,"pair_Nvertices_PLOT_abseta_bin%d_&_pt_bin%d", be, b); 
            sprintf(baff,"%s_vs_vtx_%s_pt_bin%d", idname.Data(), (be ? "endcaps":"barrel"), b); 
            single(fit_vtx,baff,buff);
        } }
        doLogX = true;
        reflbl = "Data, 1 Vtx";
        preliminary = "";
        datalbl = "Data, 2 Vtx";
        refstackNamed(fit_vtx, idname+"_vs_2vtx_pt_barrel",  "pt_PLOT_abseta_bin0_&_pair_Nvertices_bin1", "pt_PLOT_abseta_bin0_&_pair_Nvertices_bin0");
        refstackNamed(fit_vtx, idname+"_vs_2vtx_pt_endcaps", "pt_PLOT_abseta_bin1_&_pair_Nvertices_bin1", "pt_PLOT_abseta_bin1_&_pair_Nvertices_bin0");
        datalbl = "Data, 3 Vtx";
        refstackNamed(fit_vtx, idname+"_vs_3vtx_pt_barrel",  "pt_PLOT_abseta_bin0_&_pair_Nvertices_bin2", "pt_PLOT_abseta_bin0_&_pair_Nvertices_bin0");
        refstackNamed(fit_vtx, idname+"_vs_3vtx_pt_endcaps", "pt_PLOT_abseta_bin1_&_pair_Nvertices_bin2", "pt_PLOT_abseta_bin1_&_pair_Nvertices_bin0");
    }

    if (1 || ref == 0) {
        doCanvas(fit_pt, 1, 12, idname+"_barrel_pt_%d",  "abseta_bin0__pt_bin%d_");
        doCanvas(fit_pt, 1, 12, idname+"_endcaps_pt_%d", "abseta_bin1__pt_bin%d_");
    }
}



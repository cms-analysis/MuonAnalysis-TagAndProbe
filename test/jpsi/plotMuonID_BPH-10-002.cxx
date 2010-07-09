/**
  USAGE: 
    DATA-only:    root.exe -b -l -q TnP_ICHEP_MuonID_data_all.root "plotMuonID.cxx(\"data_all\")"
    MC data-like: root.exe -b -l -q TnP_ICHEP_MuonID_datalike_mc.root  "plotMuonID.cxx(\"datalike_mc\")"
    DATA+MC:      root.exe -b -l -q TnP_ICHEP_MuonID_data_all.root TnP_ICHEP_MuonID_datalike_mc.root  "plotMuonID.cxx(\"data_vs_mc\")"

    Tracker Probes vs Calo Probes:
        root.exe -b -l -q TnP_ICHEP_MuonID_FromTK_data_all.root TnP_ICHEP_MuonID_data_all.root  "plotMuonID_ICHEP.cxx(\"data_all\",2)"
        root.exe -b -l -q TnP_ICHEP_MuonID_FromTK_datalike_mc.root TnP_ICHEP_MuonID_datalike_mc.root  "plotMuonID_ICHEP.cxx(\"datalike_mc\",2)"

  REQUIRES:
   1) mkdir -p plots_bph_10_002_dev/muonid/ plots_bph_10_002_dev/muonid_tk/ plots_bph_10_002_dev/muonid_tk_vs_cal/
   2) provide a suitable "tdrStyle.cc" macro or similar
      (by default, it's taken from ~/cpp/tdrstyle.cc;
       if you need one you might want to grab ~gpetrucc/cpp/tdrstyle.cc)
*/
#include <TCanvas.h>
#include <TPad.h>
#include "plotUtil.cxx"
TString prefix = "plots_bph_10_002_dev/muonid/";
TString basedir  = "histoMuFromTk";
TString basedir2 = "histoMuFromTk";

TFile *ref = 0;

TCanvas *c1 = 0;
void plotMuonID_BPH_10_002(TString scenario, int fromTk=0) {
    prefix = prefix+scenario+"/";
    gSystem->mkdir(prefix,true);

    gROOT->ProcessLine(".x /afs/cern.ch/user/g/gpetrucc/cpp/tdrstyle.cc");
    gStyle->SetOptStat(0);
    c1 = new TCanvas("c1","c1");

    if (gROOT->GetListOfFiles()->GetEntries() == 2) {
        ref = (TFile *) gROOT->GetListOfFiles()->At(1);
        ((TFile*) gROOT->GetListOfFiles()->At(0))->cd();
    }

    datalbl = "Data, 77 nb^{-1}";
    preliminary = "";
    doRatioPlot = false; doDiffPlot = false;
    doPdf = false;
    doLogX = true;
    plotMuonIDData();
}

void plotMuonIDData() {
    char *ids[3] = { "Glb", "POG_TMLSAT", "TM" };
    for (size_t i = 0; i < 3; ++i) {
        TString idname(ids[i]);

        TDirectory *fit_pt = gFile->GetDirectory(basedir+"/"+idname+"_pt_abseta/");
        if (ref != 0) {
            TDirectory *ref_pt = ref->GetDirectory(basedir2+"/"+idname+"_pt_abseta/");
            refstack(fit_pt, ref_pt, idname+"_pt_barrel",  "pt_PLOT_abseta_bin0_");
            refstack(fit_pt, ref_pt, idname+"_pt_endcaps", "pt_PLOT_abseta_bin1_");
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

        if ((ref == 0) && (basedir2 == basedir)) {
            doCanvas(fit_pt, 1, 8, idname+"_barrel_pt_%d",   "abseta_bin0__pt_bin%d_");
            doCanvas(fit_pt, 1, 8, idname+"_endcaps_pt_%d",  "abseta_bin1__pt_bin%d_");
        }
    }

}



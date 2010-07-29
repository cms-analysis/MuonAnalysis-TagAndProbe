/**
  USAGE: 
    DATA-only:    root.exe -b -l -q TnP_ICHEP_MuonID_data_all.root "plotMuonID.cxx(\"data_all\")"
    MC data-like: root.exe -b -l -q TnP_ICHEP_MuonID_datalike_mc.root  "plotMuonID.cxx(\"datalike_mc\")"
    DATA+MC:      root.exe -b -l -q TnP_ICHEP_MuonID_data_all.root TnP_ICHEP_MuonID_datalike_mc.root  "plotMuonID.cxx(\"data_vs_mc\")"

    Tracker Probes vs Calo Probes:
        root.exe -b -l -q TnP_ICHEP_MuonID_FromTK_data_all.root TnP_ICHEP_MuonID_data_all.root  "plotMuonID_ICHEP.cxx(\"data_all\",2)"
        root.exe -b -l -q TnP_ICHEP_MuonID_FromTK_datalike_mc.root TnP_ICHEP_MuonID_datalike_mc.root  "plotMuonID_ICHEP.cxx(\"datalike_mc\",2)"

  REQUIRES:
   1) mkdir -p plots_ichep_dev/muonid/ plots_ichep_dev/muonid_tk/ plots_ichep_dev/muonid_tk_vs_cal/
   2) provide a suitable "tdrStyle.cc" macro or similar
      (by default, it's taken from ~/cpp/tdrstyle.cc;
       if you need one you might want to grab ~gpetrucc/cpp/tdrstyle.cc)
*/
#include <TCanvas.h>
#include <TPad.h>
#include "plotUtil.cxx"
TString prefix = "plots_ichep_dev/muonid/";
TString basedir  = "histoMuFromCal";
TString basedir2 = "histoMuFromCal";

TFile *ref = 0;

TCanvas *c1 = 0;
void plotMuonID_ICHEP(TString scenario, int fromTk=0) {
    if (fromTk == 1) {
        prefix = "plots_ichep_dev/muonid_tk/";
        basedir = "histoMuFromTk";
        basedir2 = "histoMuFromTk";
    } else if (fromTk == 2) {
        prefix = "plots_ichep_dev/muonid_tk_vs_cal/";
        basedir   = "histoMuFromTk";
        basedir2  = "histoMuFromCal";
        datalbl = "Trk Probes";
        reflbl = "Cal Probes";
    }

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
    datalbl = "Data, 84 nb^{-1}";
    reflbl  = "Simulation";
    preliminary = "CMS Preliminary,   #sqrt{s} = 7 TeV";
    plotMuonIDData();
}

void plotMuonIDData() {
    retitle = "Efficiency";

    char *ids[3]    = { "POG_Glb", "POG_TMLSAT", "VBTFLike" };
    char *titles[3] = { "Global",  "Soft",       "Tight" };
    for (size_t i = 0; i < 3; ++i) {
        TString idname(ids[i]);
        retitle = TString(titles[i])+" muon efficiency";
        TDirectory *fit_pt_eta = gFile->GetDirectory(basedir+"/"+idname+"_pt_abseta/");
        if (ref != 0) {
            TDirectory *ref_pt_eta = ref->GetDirectory(basedir2+"/"+idname+"_pt_abseta/");
            extraSpam = "        |#eta| < 1.2"; refstack(fit_pt_eta, ref_pt_eta, idname+"_pt_barrel",  "pt_PLOT_abseta_bin0_");
            extraSpam = "  1.2 < |#eta| < 2.4"; refstack(fit_pt_eta, ref_pt_eta, idname+"_pt_endcaps", "pt_PLOT_abseta_bin1_");
            TDirectory *mc_pt_eta  = ref->GetDirectory(basedir+"/"+idname+"_pt_abseta_mcTrue/");
            if (0 && mc_pt_eta) {
                extraSpam = "            |#eta| < 1.2"; refstack3(fit_pt_eta, ref_pt_eta, mc_pt_eta, idname+"_pt_barrel_3",  "pt_PLOT_abseta_bin0_");
                extraSpam = "       1.2 < |#eta| < 2.4"; refstack3(fit_pt_eta, ref_pt_eta, mc_pt_eta, idname+"_pt_endcaps_3", "pt_PLOT_abseta_bin1_");
            }
        } else {
            TDirectory *mc_pt_eta  = gFile->GetDirectory(basedir+"/"+idname+"_pt_abseta_mcTrue/");
            if (mc_pt_eta) {
                datalbl = "T&P fit"; reflbl = "Sim. truth";
                extraSpam = "        |#eta| < 1.2"; mcstack(fit_pt_eta, mc_pt_eta, idname+"_pt_barrel",  "pt_PLOT_abseta_bin0_");
                extraSpam = "  1.2 < |#eta| < 2.4"; mcstack(fit_pt_eta, mc_pt_eta, idname+"_pt_endcaps", "pt_PLOT_abseta_bin1_");
            } else {
                extraSpam = "        |#eta| < 1.2"; single(fit_pt_eta, idname+"_pt_barrel",  "pt_PLOT_abseta_bin0_");
                extraSpam = "  1.2 < |#eta| < 2.4"; single(fit_pt_eta, idname+"_pt_endcaps", "pt_PLOT_abseta_bin1_");
            }
        }

        TDirectory *fit_vtx = gFile->GetDirectory(basedir+"/"+idname+"_vtx/");
        if (0 && fit_vtx) {
            char buff[1255], baff[1255];
            for (int b = 0; b < 3; ++b) { for (int be = 0; be <= 1; ++be) {
                sprintf(buff,"pair_Nvertices_PLOT_abseta_bin%d_&_pt_bin%d", be, b); 
                sprintf(baff,"%s_vs_vtx_%s_pt_bin%d", idname.Data(), (be ? "endcaps":"barrel"), b); 
                single(fit_vtx,baff,buff);
            } }
        }

        if (1 || (basedir2 == basedir) && (ref == 0)) {
            doCanvas(fit_pt_eta, 1, 7, idname+"_barrel_pt_%d",   "abseta_bin0__pt_bin%d_");
            doCanvas(fit_pt_eta, 1, 7, idname+"_endcaps_pt_%d",  "abseta_bin1__pt_bin%d_");
        }
    }

}



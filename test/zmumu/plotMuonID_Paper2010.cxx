#include <TCanvas.h>
#include <TPad.h>
#include "../jpsi/plotUtil.cxx"
TString prefix = "plots_dev/muonid/";
TString basedir  = "tpTree";

TFile *ref = 0;

TCanvas *c1 = 0;
void plotMuonID_Paper2010(TString scenario="data") {

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
    doPdf = false;
    doSquare = true; yMin = 0.7; yMax = 1.048;
    datalbl = "Data, 2010B";
    reflbl  = "Simulation";
    preliminary = ""; //CMS Preliminary,   #sqrt{s} = 7 TeV";
    plotMuonIDData();
}

void plotMuonIDData() {
    retitle = "Efficiency";

    const int nids  = 5;
    char *ids[nids]    = { "Glb",    "TMOST",   "VBTF",  "VBTFold",    "PF" };
    char *titles[nids] = { "Global",  "Soft",   "Tight", "Old Tight",  "PF" };

   
    for (size_t i = 0; i < nids; ++i) {
        TString idname(ids[i]);
        retitle = TString(titles[i])+" muon efficiency";
        TDirectory *fit_pt_eta = gFile->GetDirectory(basedir+"/"+idname+"_pt_abseta/");
        if (fit_pt_eta == 0) { if (i == 0) { gFile->GetDirectory(basedir)->ls(); } ; continue; }
        if (ref != 0) {
            TDirectory *ref_pt_eta = ref->GetDirectory(basedir+"/"+idname+"_pt_abseta/");
            extraSpam = "        |#eta| < 1.2"; refstack(fit_pt_eta, ref_pt_eta, idname+"_pt_barrel",  "pt_PLOT_abseta_bin0");
            extraSpam = "  1.2 < |#eta| < 2.4"; refstack(fit_pt_eta, ref_pt_eta, idname+"_pt_endcaps", "pt_PLOT_abseta_bin1");
            TDirectory *mc_pt_eta  = ref->GetDirectory(basedir+"/"+idname+"_pt_abseta"+"_mcTrue/");
            if (0 && mc_pt_eta) {
                extraSpam = "             |#eta| < 1.2"; refstack3(fit_pt_eta, ref_pt_eta, mc_pt_eta, idname+"_pt_barrel_3",  "pt_PLOT_abseta_bin0");
                extraSpam = "       1.2 < |#eta| < 2.4"; refstack3(fit_pt_eta, ref_pt_eta, mc_pt_eta, idname+"_pt_endcaps_3", "pt_PLOT_abseta_bin1");
            }
        } else {
            TDirectory *mc_pt_eta  = gFile->GetDirectory(basedir+"/"+idname+"_pt_abseta_mcTrue/");
            if (mc_pt_eta) {
                datalbl = "T&P fit"; reflbl = "Sim. truth";
                extraSpam = "        |#eta| < 1.2"; mcstack(fit_pt_eta, mc_pt_eta, idname+"_pt_barrel",  "pt_PLOT_abseta_bin0");
                extraSpam = "  1.2 < |#eta| < 2.4"; mcstack(fit_pt_eta, mc_pt_eta, idname+"_pt_endcaps", "pt_PLOT_abseta_bin1");
            } else {
                extraSpam = "        |#eta| < 1.2"; single(fit_pt_eta, idname+"_pt_barrel",  "pt_PLOT_abseta_bin0");
                extraSpam = "  1.2 < |#eta| < 2.4"; single(fit_pt_eta, idname+"_pt_endcaps", "pt_PLOT_abseta_bin1");
            }
        }
        if (ref == 0) {
            doCanvas(fit_pt_eta, 1, 27, idname+"_barrel_pt_%d",   "abseta_bin0__pt_bin%d_");
            doCanvas(fit_pt_eta, 1, 27, idname+"_endcaps_pt_%d",  "abseta_bin1__pt_bin%d_");
        }

        TDirectory *fit_eta = gFile->GetDirectory(basedir+"/"+idname+"_eta/");
        if (ref != 0) {
            TDirectory *ref_eta = ref->GetDirectory(basedir+"/"+idname+"_eta/");
            refstack(fit_eta, ref_eta, idname+"_eta",  "eta_PLOT_");
            TDirectory *mc_eta  = ref->GetDirectory(basedir+"/"+idname+"_eta_mcTrue/");
            if (0 && mc_eta) {
                refstack3(fit_eta, ref_eta, mc_eta, idname+"_eta_3",  "eta_PLOT_");
            }
        } else {
            TDirectory *mc_eta  = gFile->GetDirectory(basedir+"/"+idname+"_eta_mcTrue/");
            if (mc_eta) {
                datalbl = "T&P fit"; reflbl = "Sim. truth";
                mcstack(fit_eta, mc_eta, idname+"_eta",  "eta_PLOT_");
            } else {
                single(fit_eta, idname+"_eta",  "eta_PLOT_");
            }
        }
        if (ref == 0) {
            doCanvas(fit_eta, 1, 27, idname+"_eta_%d",   "eta_bin%d_");
        }

        TDirectory *fit_vtx = gFile->GetDirectory(basedir+"/"+idname+"_vtx/");
        if (ref != 0) {
            TDirectory *ref_vtx = ref->GetDirectory(basedir+"/"+idname+"_vtx/");
            refstack(fit_vtx, ref_vtx, idname+"_vtx",  "tag_nVertices_PLOT_");
            TDirectory *mc_vtx  = ref->GetDirectory(basedir+"/"+idname+"_vtx_mcTrue/");
            if (0 && mc_vtx) {
                refstack3(fit_vtx, ref_vtx, mc_vtx, idname+"_vtx_3",  "tag_nVertices_PLOT_");
            }
        } else {
            TDirectory *mc_vtx  = gFile->GetDirectory(basedir+"/"+idname+"_vtx_mcTrue/");
            if (mc_vtx) {
                datalbl = "T&P fit"; reflbl = "Sim. truth";
                mcstack(fit_vtx, mc_vtx, idname+"_vtx",  "tag_nVertices_PLOT_");
            } else {
                single(fit_vtx, idname+"_vtx",  "tag_nVertices_PLOT_");
            }
        }
        if (ref == 0) {
            doCanvas(fit_vtx, 1, 27, idname+"_vtx_%d",   "tag_nVertices_bin%d_");
        }

    }
}


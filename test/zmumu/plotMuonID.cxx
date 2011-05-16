#include <TCanvas.h>
#include <TPad.h>
#include "../jpsi/plotUtil.cxx"

void plotMuonIDData() ;

void plotMuonID(TString scenario="data") {

    prefix = prefix+scenario+"/";
    gSystem->mkdir(prefix,true);

    gROOT->ProcessLine(".x /afs/cern.ch/user/g/gpetrucc/cpp/tdrstyle.cc");
    gStyle->SetOptStat(0);
    c1 = new TCanvas("c1","c1");

    if (gROOT->GetListOfFiles()->GetEntries() == 2) {
        ref = (TFile *) gROOT->GetListOfFiles()->At(1);
        ((TFile*) gROOT->GetListOfFiles()->At(0))->cd();
    }

    //fOut = TFile::Open(prefix+"plots.root", "RECREATE");
    ((TFile*) gROOT->GetListOfFiles()->At(0))->cd();

    doRatioPlot = false;
    doDiffPlot = false;
    doPdf = true;
    doSquare = true; yMin = 0.7; yMax = 1.048;
    datalbl = "Data, 2011";
    reflbl  = "Sim., 2011";
    preliminary = "CMS Preliminary,   #sqrt{s} = 7 TeV";

    if (scenario.Contains("_vs_2010")) {
        datalbl = "Data, 2011A";
        reflbl  = "Data, 2010B";
        yMinD = -0.2; yMaxD = 0.2;
        yMinR =  0.8; yMaxR = 1.2;
    }
    if (scenario.Contains("v2_vs_v1")) {
        datalbl = "2011A v2";
        reflbl  = "2011A v1";
        yMinD = -0.2; yMaxD = 0.2;
        yMinR =  0.8; yMaxR = 1.2;
    }



    plotMuonIDData();
}

void plotMuonIDData() {
    retitle = "Efficiency";

    const int nids  = 12;
    const char *ids[nids]    = { "Glb",    "TMOST",   "VBTF",  "PF", "TM",      "TMA",        "TMOSL",  "Isol_from_VBTF", "Mu24_from_VBTF", "DoubleMu7_from_VBTF", "Mu8_forEMu_from_VBTF", "Mu17_forEMu_from_VBTF" };
    const char *titles[nids] = { "Global",  "Soft",   "Tight", "PF", "Tracker", "TrackerArb", "Softer", "Isolation",      "HLT Mu24",       "HLT DoubleMu7 leg",   "HLT Mu8 leg", "HLT Mu17 leg" };

   
    for (int i = 0; i < nids; ++i) {
        TString idname(ids[i]);
        retitle = TString(titles[i])+" muon efficiency";
        if (retitle == "Isolation muon efficiency") retitle = "Isolation efficiency";
        if (retitle.Contains("HLT")) retitle = TString(titles[i])+" efficiency";
        TDirectory *fit_pt_eta = gFile->GetDirectory(basedir+"/"+idname+"_pt_abseta/");
        if (fit_pt_eta == 0) { /*if (i == 0) { gFile->GetDirectory(basedir)->ls(); } ;*/ continue; }
        yMin = 0.0; yMax = 1.1;
        if (ref != 0) {
            TDirectory *ref_pt_eta = ref->GetDirectory(basedir+"/"+idname+"_pt_abseta/");
            extraSpam = "        |#eta| < 1.2"; refstack(fit_pt_eta, ref_pt_eta, idname+"_pt_barrel",  "pt_PLOT_abseta_bin0");
            extraSpam = "  1.2 < |#eta| < 2.4"; refstack(fit_pt_eta, ref_pt_eta, idname+"_pt_endcaps", "pt_PLOT_abseta_bin1");
            TDirectory *mc_pt_eta  = ref->GetDirectory(basedir+"/"+idname+"_pt_abseta"+"_mcTrue/");
            if (mc_pt_eta) {
                extraSpam = "             |#eta| < 1.2"; refstack3(fit_pt_eta, ref_pt_eta, mc_pt_eta, idname+"_pt_barrel_3",  "pt_PLOT_abseta_bin0");
                extraSpam = "       1.2 < |#eta| < 2.4"; refstack3(fit_pt_eta, ref_pt_eta, mc_pt_eta, idname+"_pt_endcaps_3", "pt_PLOT_abseta_bin1");
            }
        } else {
            TDirectory *mc_pt_eta  = gFile->GetDirectory(basedir+"/"+idname+"_pt_abseta_mcTrue/");
            if (mc_pt_eta) {
                TString databk = datalbl, refbk = reflbl;
                datalbl = "T&P fit"; reflbl = "Sim. truth";
                extraSpam = "        |#eta| < 1.2"; mcstack(fit_pt_eta, mc_pt_eta, idname+"_pt_barrel",  "pt_PLOT_abseta_bin0");
                extraSpam = "  1.2 < |#eta| < 2.4"; mcstack(fit_pt_eta, mc_pt_eta, idname+"_pt_endcaps", "pt_PLOT_abseta_bin1");
                reflbl = refbk; datalbl = databk;
            } else {
                extraSpam = "        |#eta| < 1.2"; single(fit_pt_eta, idname+"_pt_barrel",  "pt_PLOT_abseta_bin0");
                extraSpam = "  1.2 < |#eta| < 2.4"; single(fit_pt_eta, idname+"_pt_endcaps", "pt_PLOT_abseta_bin1");
            }
        }
        //doCanvas(fit_pt_eta, 1, 7, idname+"_barrel_pt_%d",   "abseta_bin0__pt_bin%d_");
        //doCanvas(fit_pt_eta, 1, 7, idname+"_endcaps_pt_%d",  "abseta_bin1__pt_bin%d_");

        yMin = 0.7; yMax = 1.049;
        TDirectory *fit_eta = gFile->GetDirectory(basedir+"/"+idname+"_eta/");
        if (idname.Contains("Mu24")) { yMin = 0.0; yMax = 1.1; }
        if (fit_eta) {
            if (ref != 0) {
                extraSpam = "    p_{T} > 20 GeV"; retitleX = "muon #eta";
                TDirectory *ref_eta = ref->GetDirectory(basedir+"/"+idname+"_eta/");
                refstack(fit_eta, ref_eta, idname+"_eta",  "eta_PLOT_");
                TDirectory *mc_eta  = ref->GetDirectory(basedir+"/"+idname+"_eta_mcTrue/");
                if (mc_eta) {
                    refstack3(fit_eta, ref_eta, mc_eta, idname+"_eta_3",  "eta_PLOT_");
                }
            } else {
                TDirectory *mc_eta  = gFile->GetDirectory(basedir+"/"+idname+"_eta_mcTrue/");
                if (mc_eta) {
                    TString databk = datalbl, refbk = reflbl;
                    datalbl = "T&P fit"; reflbl = "Sim. truth";
                    mcstack(fit_eta, mc_eta, idname+"_eta",  "eta_PLOT_");
                    reflbl = refbk; datalbl = databk;
                } else {
                    single(fit_eta, idname+"_eta",  "eta_PLOT_");
                }
            }
            //doCanvas(fit_eta, 1, 27, idname+"_eta_%d",   "eta_bin%d_");
        }

        yMin = 0.0; yMax = 1.1;
        fit_eta = gFile->GetDirectory(basedir+"/"+idname+"_eta_fine/");
        if (fit_eta) {
            if (ref != 0) {
                extraSpam = "    p_{T} > 20 GeV"; retitleX = "muon #eta";
                TDirectory *ref_eta = ref->GetDirectory(basedir+"/"+idname+"_eta_fine/");
                refstack(fit_eta, ref_eta, idname+"_eta_fine",  "eta_PLOT_");
                TDirectory *mc_eta  = ref->GetDirectory(basedir+"/"+idname+"_eta_fine_mcTrue/");
                if (mc_eta) {
                    refstack3(fit_eta, ref_eta, mc_eta, idname+"_eta_fine_3",  "eta_PLOT_");
                }
            } else {
                TDirectory *mc_eta  = gFile->GetDirectory(basedir+"/"+idname+"_eta_fine_mcTrue/");
                if (mc_eta) {
                    TString databk = datalbl, refbk = reflbl;
                    datalbl = "T&P fit"; reflbl = "Sim. truth";
                    mcstack(fit_eta, mc_eta, idname+"_eta_fine",  "eta_PLOT_");
                    reflbl = refbk; datalbl = databk;
                } else {
                    single(fit_eta, idname+"_eta_fine",  "eta_PLOT_");
                }
            }
            if (ref == 0) {
                //doCanvas(fit_eta, 1, 27, idname+"_eta_fine_%d",   "eta_bin%d_");
            }
        }


        yMin = 0.85; yMax = 1.019;
        TDirectory *fit_vtx = gFile->GetDirectory(basedir+"/"+idname+"_vtx/");
        if (idname.Contains("Mu24")) { yMin = 0.0; yMax = 1.1; }
        if (fit_vtx) {
            TDirectory *ref_vtx = ref ? ref->GetDirectory(basedir+"/"+idname+"_vtx/") : 0;
            if (ref_vtx) {
                refstack(fit_vtx, ref_vtx, idname+"_vtx",  "tag_nVertices_PLOT_");
                TDirectory *mc_vtx  = ref->GetDirectory(basedir+"/"+idname+"_vtx_mcTrue/");
                if (mc_vtx) {
                    refstack3(fit_vtx, ref_vtx, mc_vtx, idname+"_vtx_3",  "tag_nVertices_PLOT_");
                }
            } else {
                TDirectory *mc_vtx  = gFile->GetDirectory(basedir+"/"+idname+"_vtx_mcTrue/");
                if (mc_vtx) {
                    TString databk = datalbl, refbk = reflbl;
                    datalbl = "T&P fit"; reflbl = "Sim. truth";
                    mcstack(fit_vtx, mc_vtx, idname+"_vtx",  "tag_nVertices_PLOT_");
                    reflbl = refbk; datalbl = databk;
                } else {
                    single(fit_vtx, idname+"_vtx",  "tag_nVertices_PLOT_");
                }
            }
            if (ref == 0) doCanvas(fit_vtx, 1, 27, idname+"_vtx_%d",   "tag_nVertices_bin%d_");
        }

        yMin = 0.85; yMax = 1.019;
        fit_vtx = gFile->GetDirectory(basedir+"/"+idname+"_vtxDA/");
        if (fit_vtx) {
            retitleX = "number of vertices (DA 100#mum)";
            TDirectory *ref_vtx = ref ? ref->GetDirectory(basedir+"/"+idname+"_vtxDA/") : 0;
            if (ref_vtx) {
                refstack(fit_vtx, ref_vtx, idname+"_vtxDA",  "tag_nVerticesDA_PLOT_");
                TDirectory *mc_vtx  = ref->GetDirectory(basedir+"/"+idname+"_vtxDA_mcTrue/");
                if (mc_vtx) {
                    refstack3(fit_vtx, ref_vtx, mc_vtx, idname+"_vtxDA_3",  "tag_nVerticesDA_PLOT_");
                }
            } else {
                TDirectory *mc_vtx  = gFile->GetDirectory(basedir+"/"+idname+"_vtxDA_mcTrue/");
                if (mc_vtx) {
                    TString databk = datalbl, refbk = reflbl;
                    datalbl = "T&P fit"; reflbl = "Sim. truth";
                    mcstack(fit_vtx, mc_vtx, idname+"_vtxDA",  "tag_nVerticesDA_PLOT_");
                    reflbl = refbk; datalbl = databk;
                } else {
                    single(fit_vtx, idname+"_vtxDA",  "tag_nVerticesDA_PLOT_");
                }
            }
            if (ref == 0) doCanvas(fit_vtx, 1, 27, idname+"_vtxDA_%d",   "tag_nVerticesDA_bin%d_");
        }


        yMin = 0.85; yMax = 1.019;
        TDirectory *fit_charge = gFile->GetDirectory(basedir+"/"+idname+"_charge/");
        if (fit_charge) {
            TDirectory *ref_charge = ref ? ref->GetDirectory(basedir+"/"+idname+"_charge/") : 0;
            if (ref_charge != 0) {
                refstack(fit_charge, ref_charge, idname+"_charge",  "charge_PLOT_");
                TDirectory *mc_charge  = ref->GetDirectory(basedir+"/"+idname+"_charge_mcTrue/");
                if (mc_charge) {
                    refstack3(fit_charge, ref_charge, mc_charge, idname+"_charge_3",  "charge_PLOT_");
                }
            } else {
                TDirectory *mc_charge  = gFile->GetDirectory(basedir+"/"+idname+"_charge_mcTrue/");
                if (mc_charge) {
                    TString databk = datalbl, refbk = reflbl;
                    datalbl = "T&P fit"; reflbl = "Sim. truth";
                    mcstack(fit_charge, mc_charge, idname+"_charge",  "charge_PLOT_");
                    reflbl = refbk; datalbl = databk;
                } else {
                    single(fit_charge, idname+"_charge",  "charge_PLOT_");
                }
            }
            if (ref == 0) {
                //doCanvas(fit_charge, 1, 2, idname+"_charge_%d",   "charge_bin%d_");
            }
        }

        yMin = 0.941; yMax = 1.009;
        if (idname.Contains("VBTF")) { yMin -= 0.02; yMax  -= 0.02; }

        TDirectory *fit_overall = gFile->GetDirectory(basedir+"/"+idname+"_overall/");
        if (fit_overall) {
            TDirectory *ref_overall = ref ? ref->GetDirectory(basedir+"/"+idname+"_overall/") : 0;
            if (ref_overall != 0) {
                refstack(fit_overall, ref_overall, idname+"_overall",  "abseta_PLOT_");
                TDirectory *mc_overall  = ref->GetDirectory(basedir+"/"+idname+"_overall_mcTrue/");
                if (mc_overall) {
                    refstack3(fit_overall, ref_overall, mc_overall, idname+"_overall_3",  "abseta_PLOT_");
                }
            } else {
                TDirectory *mc_overall  = gFile->GetDirectory(basedir+"/"+idname+"_overall_mcTrue/");
                if (mc_overall) {
                    TString databk = datalbl, refbk = reflbl;
                    datalbl = "T&P fit"; reflbl = "Sim. truth";
                    mcstack(fit_overall, mc_overall, idname+"_overall",  "abseta_PLOT_");
                    reflbl = refbk; datalbl = databk;
                } else {
                    single(fit_overall, idname+"_overall",  "abseta_PLOT_");
                }
            }
            //doCanvas(fit_overall, 1, 1, idname+"_overall_%d",   "abseta_bin%d_");
        }

        yMin = 0.85; yMax = 1.019;
        TDirectory *fit_overall_abseta = gFile->GetDirectory(basedir+"/"+idname+"_overall_abseta/");
        if (fit_overall_abseta) {
            TDirectory *ref_overall_abseta = ref ? ref->GetDirectory(basedir+"/"+idname+"_overall_abseta/") : 0;
            if (ref_overall_abseta) {
                extraSpam = "    p_{T} > 20 GeV"; retitleX = "muon |#eta|";
                refstack(fit_overall_abseta, ref_overall_abseta, idname+"_overall_abseta",  "abseta_PLOT_");
            } else {
                TDirectory *mc_overall_abseta  = gFile->GetDirectory(basedir+"/"+idname+"_overall_abseta_mcTrue/");
                if (mc_overall_abseta) {
                    TString databk = datalbl, refbk = reflbl;
                    datalbl = "T&P fit"; reflbl = "Sim. truth";
                    mcstack(fit_overall_abseta, mc_overall_abseta, idname+"_overall_abseta",  "abseta_PLOT_");
                    reflbl = refbk; datalbl = databk;
                } else {
                    single(fit_overall_abseta, idname+"_overall_abseta",  "abseta_PLOT_");
                }
            }
            doCanvas(fit_overall_abseta, 1, 27, idname+"_overall_abseta_%d",   "abseta_bin%d_");
        }
        TDirectory *fit_overall_endcaps21 = gFile->GetDirectory(basedir+"/"+idname+"_overall_endcaps21/");
        if (fit_overall_endcaps21) {
            TDirectory *ref_overall_endcaps21 = ref ? ref->GetDirectory(basedir+"/"+idname+"_overall_endcaps21/") : 0;
            if (ref_overall_endcaps21) {
                extraSpam = "    p_{T} > 20 GeV"; retitleX = "muon |#eta|";
                refstack(fit_overall_endcaps21, ref_overall_endcaps21, idname+"_overall_endcaps21",  "abseta_PLOT_");
            } else {
                single(fit_overall_endcaps21, idname+"_overall_endcaps21",  "abseta_PLOT_");
            }
            doCanvas(fit_overall_endcaps21, 1, 27, idname+"_overall_endcaps21_%d",   "abseta_bin%d_");
        }
    }
}


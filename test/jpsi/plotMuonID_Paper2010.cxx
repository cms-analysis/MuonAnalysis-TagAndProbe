#include <TCanvas.h>
#include <TPad.h>
#include "plotUtil.cxx"
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

    fOut = TFile::Open(prefix+"plots.root", "RECREATE");
    ((TFile*) gROOT->GetListOfFiles()->At(0))->cd();

    doFillMC = true;
    doRatioPlot = true;
    doDiffPlot  = true;
    doPdf = true;
    doSquare = true; yMin = 0; yMax = 1.1;
    datalbl = "Data, 2010";
    reflbl  = "Simulation";
    if (scenario.Contains("tk_vs_calo")) {
        datalbl = "All tracks";
        reflbl  = "MIP tracks";
        yMinD = -0.2; yMaxD = 0.2;
        yMinR =  0.8; yMaxR = 1.2;
    }
    if (scenario.Contains("39X_vs_38X")) {
        TString what = (scenario.Contains("mc") ? " MC" : " Data");
        datalbl = "3.9.X"+what;
        reflbl  = "3.8.X"+what;
        yMinD = -0.2; yMaxD = 0.2;
        yMinR =  0.8; yMaxR = 1.2;
    }

    if (scenario.Contains("promopt_vs_bdecay") || scenario.Contains("prompt_vs_bdecay") ) {
        datalbl = "B decays";
        reflbl  = "Prompt";
        yMinD = -0.1; yMaxD = 0.1;
        yMinR =  0.9; yMaxR = 1.1;
    }
    if (scenario.Contains("cbq_vs_gaussexp")) {
        if (TString(gFile->GetName()).Contains("gaussExpo")) {
            datalbl  = "Gaus + Exp";
            reflbl = "CB + Pol2";
        } else {
            datalbl = "CB + Pol2";
            reflbl  = "Gaus + Exp";
        }
        yMinD = -0.2; yMaxD = 0.2;
        yMinR =  0.8; yMaxR = 1.2;
    }

    if (ref == 0 && scenario.Contains("mc")) {
        yMinD = -0.04; yMaxD = 0.04;
        yMinR =  0.96; yMaxR = 1.04;
    }
    preliminary = "CMS Preliminary,   #sqrt{s} = 7 TeV";
    //doLogX = true;

    plotMuonIDData();
}

TGraphAsymmErrors *merge(TGraphAsymmErrors **fits, int ntrig, int neta, int ieta) {
    int founds[9];
    int nfound = 0;
    for (size_t j = 0; j < ntrig; ++j) {
        if (fits[j*neta + ieta] != 0 && fits[j*neta + ieta]->GetN() != 0) founds[nfound++] = j;
    }
    if (nfound == 0) return 0;
    if (nfound > 1) {
        TGraphAsymmErrors *g1 = fits[founds[0]*neta+ieta], *g2 = fits[founds[1]*neta+ieta];
        TGraphAsymmErrors *g3 = nfound > 2 ? fits[founds[2]*neta+ieta] : NULL;
        TGraphAsymmErrors *g4 = nfound > 3 ? fits[founds[3]*neta+ieta] : NULL;
        TGraphAsymmErrors *merged = mergeGraphs(g1,g2,g3,g4,Last);
        merged->GetXaxis()->SetTitle(g1->GetXaxis()->GetTitle());
        merged->GetYaxis()->SetTitle(g1->GetYaxis()->GetTitle());
        return merged;
    } else {
        return fits[founds[0]*neta+ieta];
    }

}

void plotMuonIDData() {
    retitle = "Efficiency";

    const int nids  = 4;
    char *ids[nids]    = { "TMOST",   "VBTF", "PF", "Glb"    };
    char *titles[nids] = { "Soft",   "Tight", "PF", "Global" };

    const int ntrig = 4;
    char *trig[ntrig] = { "Mu3_Track0", "Mu5_Track0", "Mu3_Track3", "Mu3_Track5" };

    const int neta = 2;
    TGraphAsymmErrors *fits[nids][ntrig][neta];
    TGraphAsymmErrors *refs[nids][ntrig][neta];
    TGraphAsymmErrors *mcts[nids][ntrig][neta];
    for (size_t i = 0; i < nids; ++i) { for (size_t j = 0; j < ntrig; ++j) { for (size_t k = 0; k < neta; ++k) { fits[i][j][k] = 0; } } }
    for (size_t i = 0; i < nids; ++i) { for (size_t j = 0; j < ntrig; ++j) { for (size_t k = 0; k < neta; ++k) { refs[i][j][k] = 0; } } }
    for (size_t i = 0; i < nids; ++i) { for (size_t j = 0; j < ntrig; ++j) { for (size_t k = 0; k < neta; ++k) { mcts[i][j][k] = 0; } } }
    
    for (size_t i = 0; i < nids; ++i) {
        for (size_t j = 0; j < ntrig; ++j) {
            yMin = 0.0; yMax = 1.1;
            TString idname(ids[i]);
            TString tname(trig[j]);
            retitle = TString(titles[i])+" muon efficiency";
            TDirectory *fit_pt_eta = gFile->GetDirectory(basedir+"/"+idname+"_pt_abseta_"+tname+"/");
            TDirectory *fit_p_eta  = gFile->GetDirectory(basedir+"/"+idname+"_p_abseta_"+tname+"/");
            if (fit_pt_eta == 0) { if (i == 0) { gFile->GetDirectory(basedir)->ls(); } ; continue; }
            TDirectory *ref_pt_eta = ref ? ref->GetDirectory(basedir+"/"+idname+"_pt_abseta_"+tname+"/") : 0;
            TDirectory *ref_p_eta = ref ? ref->GetDirectory(basedir+"/"+idname+"_p_abseta_"+tname+"/") : 0;
            if (ref_pt_eta) {
                fits[i][j][0] = getFit(fit_pt_eta, "pt_PLOT_abseta_bin0_");
                fits[i][j][1] = getFit(fit_pt_eta, "pt_PLOT_abseta_bin1_");
                refs[i][j][0] = getFit(ref_pt_eta, "pt_PLOT_abseta_bin0_");
                refs[i][j][1] = getFit(ref_pt_eta, "pt_PLOT_abseta_bin1_");
                if (neta >= 3) fits[i][j][2] = getFit(fit_p_eta, "p_PLOT_abseta_bin0_");
                if (neta >= 3) refs[i][j][2] = getFit(ref_p_eta, "p_PLOT_abseta_bin0_");
                doLogX = true;
                extraSpam = "        |#eta| < 1.2"; refstack(fit_pt_eta, ref_pt_eta, idname+"_pt_barrel_"+tname,  "pt_PLOT_abseta_bin0_");
                doLogX = false;
                //extraSpam = "  0.9 < |#eta| < 1.2"; refstack(fit_pt_eta, ref_pt_eta, idname+"_pt_overlap_"+tname, "pt_PLOT_abseta_bin1_");
                extraSpam = "  1.2 < |#eta| < 2.4"; refstack(fit_pt_eta, ref_pt_eta, idname+"_pt_endcaps_"+tname, "pt_PLOT_abseta_bin1_");
                if (neta >= 3 && fits[i][j][2] && refs[i][j][2]) {
                    extraSpam = "  1.2 < |#eta| < 2.4"; refstack(fit_p_eta, ref_p_eta, idname+"_p_endcaps_"+tname, "p_PLOT_abseta_bin0_");
                }
                TDirectory *mc_pt_eta  = ref->GetDirectory(basedir+"/"+idname+"_pt_abseta_"+tname+"_mcTrue/");
                if (mc_pt_eta) {
                    mcts[i][j][0] = getFit(mc_pt_eta, "pt_PLOT_abseta_bin0_", "cnt");
                    mcts[i][j][1] = getFit(mc_pt_eta, "pt_PLOT_abseta_bin1_", "cnt");
                }
                if (0 && mc_pt_eta) {
                    extraSpam = "             |#eta| < 1.2"; refstack3(fit_pt_eta, ref_pt_eta, mc_pt_eta, idname+"_pt_barrel_3_"+tname,  "pt_PLOT_abseta_bin0_");
                    //extraSpam = "       0.9 < |#eta| < 1.2"; refstack3(fit_pt_eta, ref_pt_eta, mc_pt_eta, idname+"_pt_overlap_3_"+tname, "pt_PLOT_abseta_bin1_");
                    extraSpam = "       1.2 < |#eta| < 2.4"; refstack3(fit_pt_eta, ref_pt_eta, mc_pt_eta, idname+"_pt_endcaps_3_"+tname, "pt_PLOT_abseta_bin1_");
                }
            } else {
                TDirectory *mc_pt_eta  = gFile->GetDirectory(basedir+"/"+idname+"_pt_abseta_"+tname+"_mcTrue/");
                if (mc_pt_eta) {
                    TString databk = datalbl, refbk = reflbl;
                    datalbl = "T&P fit"; reflbl = "Sim. truth";
                    extraSpam = "        |#eta| < 1.2"; mcstack(fit_pt_eta, mc_pt_eta, idname+"_pt_barrel_"+tname,  "pt_PLOT_abseta_bin0_");
                    extraSpam = "  1.2 < |#eta| < 2.4"; mcstack(fit_pt_eta, mc_pt_eta, idname+"_pt_endcaps_"+tname, "pt_PLOT_abseta_bin1_");
                    reflbl = refbk; datalbl = databk;
                } else {
                    extraSpam = "        |#eta| < 1.2"; fits[i][j][0] = single(fit_pt_eta, idname+"_pt_barrel_"+tname,  "pt_PLOT_abseta_bin0_");
                    extraSpam = "  1.2 < |#eta| < 2.4"; fits[i][j][1] = single(fit_pt_eta, idname+"_pt_endcaps_"+tname, "pt_PLOT_abseta_bin1_");
                    if (neta >= 3) {
                        extraSpam = "  1.2 < |#eta| < 2.4"; fits[i][j][2] = single(fit_p_eta,  idname+"_p_endcaps_"+tname,  "p_PLOT_abseta_bin1_");
                    }
                }
            }
            if (tname == "Mu3_Track3") {
                //yMin = 0.59; yMax = 1.049;
                yMin = 0.0; yMax = 1.1;
                TDirectory *fit_vtx_barrel  = gFile->GetDirectory(basedir+"/"+idname+"_vtx_barrel_" +tname+"/");
                TDirectory *fit_vtx_endcaps = gFile->GetDirectory(basedir+"/"+idname+"_vtx_endcaps_"+tname+"/");
                if (fit_vtx_barrel && fit_vtx_endcaps) {
                    if (ref) {
                        TDirectory *ref_vtx_barrel  = ref->GetDirectory(basedir+"/"+idname+"_vtx_barrel_" +tname+"/");
                        TDirectory *ref_vtx_endcaps = ref->GetDirectory(basedir+"/"+idname+"_vtx_endcaps_"+tname+"/");
                        extraSpam = "        |#eta| < 0.9"; refstack(fit_vtx_barrel,  ref_vtx_barrel,  idname+"_vtx_barrel_"+tname,  "tag_nVertices_PLOT_");
                        extraSpam = "  1.2 < |#eta| < 2.4"; refstack(fit_vtx_endcaps, ref_vtx_endcaps, idname+"_vtx_endcaps_"+tname, "tag_nVertices_PLOT_");
                    } else {
                        single(fit_vtx_barrel,  idname+"_vtx_barrel_"+tname,  "tag_nVertices_PLOT_");
                        single(fit_vtx_endcaps, idname+"_vtx_endcaps_"+tname, "tag_nVertices_PLOT_");
                        doCanvas(fit_vtx_barrel,  1, 10, idname+"_vtx_barrel_" +tname+"_vtx_%d",  ".*tag_nVertices_bin%d_");
                        doCanvas(fit_vtx_endcaps, 1, 10, idname+"_vtx_endcaps_"+tname+"_vtx_%d",  ".*tag_nVertices_bin%d_");
                    }
                }
            }

            yMin = 0.0; yMax = 1.1;
            TDirectory *fit_pt6_eta  = gFile->GetDirectory(basedir+"/"+idname+"_pt6_eta_" +tname+"/");
            if (fit_pt6_eta) {
                if (ref) {
                    TDirectory *ref_pt6_eta  = ref->GetDirectory(basedir+"/"+idname+"_pt6_eta_" +tname+"/");
                    extraSpam = "  p_{T} > 5 GeV/c"; refstack(fit_pt6_eta,  ref_pt6_eta,  idname+"_pt6_eta_"+tname,  "eta_PLOT_");
                } else {
                    single(fit_pt6_eta,  idname+"_pt6_eta_"+tname,  "eta_PLOT_");
                    doCanvas(fit_pt6_eta,  1, 20, idname+"_pt6_eta_" +tname+"_eta_%d",  ".*eta_bin%d_");
                }
            }

            yMin = 0.85; yMax = 1.019;
            TDirectory *fit_plateau_barrel  = gFile->GetDirectory(basedir+"/"+idname+"_plateau_barrel_" +tname+"/");
            TDirectory *fit_plateau_endcaps  = gFile->GetDirectory(basedir+"/"+idname+"_plateau_endcaps_" +tname+"/");
            TDirectory *fit_plateau_endcaps21  = gFile->GetDirectory(basedir+"/"+idname+"_plateau_endcaps21_" +tname+"/");
            if (fit_plateau_barrel) {
                TString tmcname = tname;
                if (tname == "Mu3_Track5" && ref != 0 && 
                        ref->GetDirectory(basedir+"/"+idname+"_plateau_barrel_"+tmcname+"/") == 0) {
                    tmcname = "Mu3_Track3";
                }
                TDirectory *ref_plateau_barrel  = ref ? ref->GetDirectory(basedir+"/"+idname+"_plateau_barrel_"  +tmcname+"/") : 0;
                TDirectory *ref_plateau_endcaps = ref ? ref->GetDirectory(basedir+"/"+idname+"_plateau_endcaps_" +tmcname+"/") : 0;
                TDirectory *ref_plateau_endcaps21 = ref ? ref->GetDirectory(basedir+"/"+idname+"_plateau_endcaps21_" +tmcname+"/") : 0;
                if (ref_plateau_barrel) {
                    extraSpam = (idname == "VBTF" ? "  p_{T} > 10 GeV/c" : "  p_{T} > 6 GeV/c"); 
                    refstack(fit_plateau_barrel,  ref_plateau_barrel,  idname+"_plateau_barrel_"+tname,  "abseta_PLOT_");
                    extraSpam = (idname == "VBTF" ? "  p_{T} > 10 GeV/c" : "  p_{T} > 4 GeV/c"); 
                    refstack(fit_plateau_endcaps,  ref_plateau_endcaps,  idname+"_plateau_endcaps_"+tname,  "abseta_PLOT_");
                } else {
                    single(fit_plateau_barrel,  idname+"_plateau_barrel_"+tname,  "abseta_PLOT_");
                    single(fit_plateau_endcaps,  idname+"_plateau_endcaps_"+tname,  "abseta_PLOT_");
                    doCanvas(fit_plateau_barrel,  1, 1, idname+"_plateau_barrel_" +tname,  "abseta_bin1_");
                    doCanvas(fit_plateau_endcaps,  1, 1, idname+"_plateau_endcaps_" +tname,  "abseta_bin1_");
                }
                if (ref_plateau_endcaps21) {
                    extraSpam = (idname == "VBTF" ? "  p_{T} > 10 GeV/c" : "  p_{T} > 4 GeV/c"); 
                    refstack(fit_plateau_endcaps21,  ref_plateau_endcaps21,  idname+"_plateau_endcaps21_"+tname,  "abseta_PLOT_");
                } else if (fit_plateau_endcaps21) {
                    single(fit_plateau_endcaps21,  idname+"_plateau_endcaps21_"+tname,  "abseta_PLOT_");
                    doCanvas(fit_plateau_endcaps21,  1, 1, idname+"_plateau_endcaps21_" +tname,  "abseta_bin1_");
                }
            }


            yMin = 0.0; yMax = 1.1;
            TDirectory *fit_pt_barrelStrict  = gFile->GetDirectory(basedir+"/"+idname+"_pt_barrel_" +tname+"/");
            if (fit_pt_barrelStrict) {
                doLogX = true;
                if (ref) {
                    TDirectory *ref_pt_barrelStrict  = ref->GetDirectory(basedir+"/"+idname+"_pt_barrel_" +tname+"/");
                    extraSpam = "        |#eta| < 0.8"; 
                    refstack(fit_pt_barrelStrict,  ref_pt_barrelStrict,  idname+"_pt_barrelStrict_"+tname,  "pt_PLOT_");
                } else {
                    single(fit_pt_barrelStrict,  idname+"_pt_barrelStrict_"+tname,  "eta_PLOT_");
                    //doCanvas(fit_pt_barrel,  1, 20, idname+"_pt_barrelStrict_" +tname+"_eta_%d",  ".*eta_bin%d_");
                }
                doLogX = false;
            }

        }
    }

    bool doFillMCSave = doFillMC;
    doFillMC = false;
    c1->cd(); c1->Clear(); 
    if (doSquare) squareCanvas(c1);
    double xmax = 0, pmax = 0; TH1F *frame = 0, *pframe = 0;
    for (size_t i = 0; i < nids; ++i) {
        int nfound = 0, jfound = -1;
        for (size_t j = 0; j < ntrig; ++j) {
            if (fits[i][j][0] && fits[i][j][1] && fits[i][j][0]->GetN() && fits[i][j][1]->GetN()) {
                xmax = TMath::Max(xmax, xmaxGraph(fits[i][j][0])); 
                xmax = TMath::Max(xmax, xmaxGraph(fits[i][j][1])); 
                if (neta >= 3) pmax = TMath::Max(pmax, xmaxGraph(fits[i][j][2])); 
                jfound = j;
            }
        }
        if (jfound == -1) continue;
        if (frame == 0) {
            frame = new TH1F("frame","frame",1,0.,xmax); 
            frame->GetYaxis()->SetRangeUser(yMin,yMax);
            frame->GetYaxis()->SetTitle(fits[i][jfound][0]->GetYaxis()->GetTitle());
            frame->GetXaxis()->SetTitle(fits[i][jfound][0]->GetXaxis()->GetTitle());
        }
        if (neta >= 3 && pframe == 0) {
            pframe = new TH1F("pframe","pframe",1,0.,pmax); 
            pframe->GetYaxis()->SetRangeUser(yMin,yMax);
            pframe->GetYaxis()->SetTitle(fits[i][jfound][2]->GetYaxis()->GetTitle());
            pframe->GetXaxis()->SetTitle(fits[i][jfound][2]->GetXaxis()->GetTitle());
        }
        int founds[ntrig];
        int colors[6] = { 1, 2, 4, 209, 6, 9 };
        for (int be = 0; be < neta; ++be) {
            frame->Draw();
            extraSpam = (be == 0 ? "        |#eta| < 1.2" : "  1.2 < |#eta| < 2.4");
            nfound = 0;
            for (size_t j = 0; j < ntrig; ++j) {
                if (fits[i][j][be]  && fits[i][j][be]->GetN()) { 
                    fits[i][j][be]->SetLineColor(colors[nfound]);
                    fits[i][j][be]->SetMarkerColor(colors[nfound]);
                    fits[i][j][be]->Draw("P");
                    founds[nfound++] = j;
                }
            }
            if (nfound == 2) doLegend(fits[i][founds[0]][be], fits[i][founds[1]][be], trig[founds[0]], trig[founds[1]]);
            else if (nfound == 3) doLegend(fits[i][founds[0]][be], fits[i][founds[1]][be], fits[i][founds[2]][be], trig[founds[0]], trig[founds[1]], trig[founds[2]]);
            else if (nfound == 4) doLegend(fits[i][founds[0]][be], fits[i][founds[1]][be], fits[i][founds[2]][be], fits[i][founds[3]][be], trig[founds[0]], trig[founds[1]], trig[founds[2]], trig[founds[3]]);
            else if (nfound > 1) { std::cerr << "ERROR: can have at most 4 triggers, sorry." << std::endl; }
            maybeLogX(c1, fits[i][founds[0]][be]); 
            TString label = TString::Format("%s_%s", ids[i], (be == 0 ? "barrel" : (be == 1 ? "endcaps" : "endcaps_p") ));
            c1->Print(prefix+label+"_all.png");
            if (doPdf) c1->Print(prefix+label+"_all.pdf");
            if (nfound > 1) {
                TGraphAsymmErrors *g1 = fits[i][founds[0]][be], *g2 = fits[i][founds[1]][be];
                TGraphAsymmErrors *g3 = nfound > 2 ? fits[i][founds[2]][be] : NULL;
                TGraphAsymmErrors *g4 = nfound > 3 ? fits[i][founds[3]][be] : NULL;
                TGraphAsymmErrors *merged = mergeGraphs(g1,g2,g3,g4);
                frame->Draw();
                merged->Draw("P");
                c1->Print(prefix+label+"_merge.png");
                if (doPdf) c1->Print(prefix+label+"_merge.pdf");
                if (fOut) fOut->WriteTObject(merged, "fit_"+label+"_merge");
            }
        }
    }
    doFillMC = doFillMCSave;
    if (ref) {
        for (size_t i = 0; i < nids; ++i) {
            for (int be = 0; be < neta; ++be) {
                TString label = TString::Format("%s_%s", ids[i], (be == 0 ? "barrel" : (be == 1 ? "endcaps" : "endcaps_p") ));
                TGraphAsymmErrors *mref = merge(refs[i], ntrig, neta, be);
                TGraphAsymmErrors *mmct = merge(mcts[i], ntrig, neta, be);
                TGraphAsymmErrors *mfit = merge(fits[i], ntrig, neta, be);
                if (mref == 0) std::cerr << "No mref for " << label << std::endl;
                if (mfit == 0) std::cerr << "No mfit for " << label << std::endl;
                if (mref ==0 || mfit == 0) continue;
                if (doFillMC) {
                    mref->SetLineColor(2);
                    mref->SetFillColor(208);
                    mref->SetLineStyle(0);
                    mref->SetMarkerColor(2);
                    mref->SetMarkerStyle(1);
                    mref->SetMarkerSize(0.4);
                } else {
                    mref->SetLineWidth(2);
                    mref->SetLineColor(kRed);
                    mref->SetMarkerColor(kRed);
                    mref->SetMarkerStyle(25);
                    mref->SetMarkerSize(2.0);
                }
                mfit->SetLineWidth(2);
                mfit->SetLineColor(kBlack);
                mfit->SetMarkerColor(kBlack);
                mfit->SetMarkerStyle(20);
                mfit->SetMarkerSize(1.6);
                if (retitle != "") 
                frame->GetYaxis()->SetTitle(TString(titles[i])+" muon efficiency");
                frame->GetXaxis()->SetTitle(be <= 1 ? "muon p_{T}  (GeV/c)" : "muon p  (GeV/c)");
                extraSpam = (be == 0 ? "        |#eta| < 1.2" : "  1.2 < |#eta| < 2.4");
                frame->Draw();
                mref->Draw(doFillMC ? "E2 SAME" : "P SAME");
                mfit->Draw("P SAME");
                doLegend(mfit, mref, datalbl, reflbl);
                maybeLogX(c1, mfit); 
                c1->Print(prefix+label+"_merge.png");
                if (doPdf) c1->Print(prefix+label+"_merge.pdf");
                if (fOut) fOut->WriteTObject(mfit, "fit_"+label+"_merge");
                if (fOut) fOut->WriteTObject(mref, "ref_"+label+"_merge");
                if (fOut && mmct) fOut->WriteTObject(mmct, "mct_"+label+"_merge");
            }
        }
    }
}


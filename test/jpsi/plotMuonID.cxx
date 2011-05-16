#include <TCanvas.h>
#include <TPad.h>
//TCanvas *c1 = 0;

#include "plotUtil.cxx"

bool doCanvases = false;

TGraphAsymmErrors *merge(TGraphAsymmErrors **fits, int ntrig, int neta, int ieta) ;
void plotMuonIDData() ;

void plotMuonID(TString scenario="data") {
    basedir  = "tpTree";
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

    doFillMC = false;
    doRatioPlot = false;
    doDiffPlot  = false;
    doPdf = true;
    doSquare = true; yMin = 0; yMax = 1.1;
    datalbl = "Data, 2011";
    reflbl  = "Sim., 2010";
    if (scenario.Contains("tk_vs_calo")) {
        datalbl = "All tracks";
        reflbl  = "MIP tracks";
        yMinD = -0.2; yMaxD = 0.2;
        yMinR =  0.8; yMaxR = 1.2;
    }
    if (scenario.Contains("39X_vs_38X")) {
        TString what = (scenario.Contains("mc") ? " Sim." : " Data");
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

    doCanvases = false; //if (scenario.Contains("data_vs_mc")) doCanvases = true;

    plotMuonIDData();
}

TGraphAsymmErrors *merge(TGraphAsymmErrors **fits, int ntrig, int neta, int ieta) {
    int founds[9];
    int nfound = 0;
    for (int j = 0; j < ntrig; ++j) {
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
    const char *ids[nids]    = { "TMOST",   "VBTF", "PF", "Glb"    };
    const char *titles[nids] = { "Soft",   "Tight", "PF", "Global" };

    const int ntrig = 5;
    const char *trig[ntrig] = { "Mu5_Track2", "Mu7_Track7", "Mu5_Track0", "Mu3_Track3", "Mu3_Track5" };

    const int neta = 2;
    TGraphAsymmErrors *fits[nids][ntrig][neta], *fits_wide[nids][ntrig][neta];
    TGraphAsymmErrors *refs[nids][ntrig][neta], *refs_wide[nids][ntrig][neta];
    TGraphAsymmErrors *mcts[nids][ntrig][neta], *mcts_wide[nids][ntrig][neta];
    for (int i = 0; i < nids; ++i) { for (int j = 0; j < ntrig; ++j) { for (int k = 0; k < neta; ++k) { fits[i][j][k] = 0; fits_wide[i][j][k] = 0; } } }
    for (int i = 0; i < nids; ++i) { for (int j = 0; j < ntrig; ++j) { for (int k = 0; k < neta; ++k) { refs[i][j][k] = 0; refs_wide[i][j][k] = 0; } } }
    for (int i = 0; i < nids; ++i) { for (int j = 0; j < ntrig; ++j) { for (int k = 0; k < neta; ++k) { mcts[i][j][k] = 0; mcts_wide[i][j][k] = 0; } } }
    
    for (int i = 0; i < nids; ++i) {
        TString idname(ids[i]);

        for (int j = 0; j < ntrig; ++j) {
            TString tname(trig[j]);
            retitle = TString(titles[i])+" muon efficiency";

            yMin = 0.0; yMax = 1.1;
            TDirectory *fit_pt_eta = gFile->GetDirectory(basedir+"/"+idname+"_pt_abseta_"+tname+"/");
            TDirectory *ref_pt_eta = ref ? ref->GetDirectory(basedir+"/"+idname+"_pt_abseta_"+tname+"/") : 0;
            if (fit_pt_eta && ref_pt_eta) {
                fits[i][j][0] = getFit(fit_pt_eta, "pt_PLOT_abseta_bin0_");
                fits[i][j][1] = getFit(fit_pt_eta, "pt_PLOT_abseta_bin1_");
                refs[i][j][0] = getFit(ref_pt_eta, "pt_PLOT_abseta_bin0_");
                refs[i][j][1] = getFit(ref_pt_eta, "pt_PLOT_abseta_bin1_");
                doLogX = true;
                extraSpam = "        |#eta| < 1.2"; refstack(fit_pt_eta, ref_pt_eta, idname+"_pt_barrel_"+tname,  "pt_PLOT_abseta_bin0_");
                doLogX = false;
                //extraSpam = "  0.9 < |#eta| < 1.2"; refstack(fit_pt_eta, ref_pt_eta, idname+"_pt_overlap_"+tname, "pt_PLOT_abseta_bin1_");
                extraSpam = "  1.2 < |#eta| < 2.4"; refstack(fit_pt_eta, ref_pt_eta, idname+"_pt_endcaps_"+tname, "pt_PLOT_abseta_bin1_");
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
            } else if (fit_pt_eta) {
                TDirectory *mc_pt_eta  = gFile->GetDirectory(basedir+"/"+idname+"_pt_abseta_"+tname+"_mcTrue/");
                if (mc_pt_eta) {
                    fits[i][j][0] = getFit(fit_pt_eta, "pt_PLOT_abseta_bin0_");
                    fits[i][j][1] = getFit(fit_pt_eta, "pt_PLOT_abseta_bin1_");
                    TString databk = datalbl, refbk = reflbl;
                    datalbl = "T&P fit"; reflbl = "Sim. truth";
                    extraSpam = "        |#eta| < 1.2"; mcstack(fit_pt_eta, mc_pt_eta, idname+"_pt_barrel_"+tname,  "pt_PLOT_abseta_bin0_");
                    extraSpam = "  1.2 < |#eta| < 2.4"; mcstack(fit_pt_eta, mc_pt_eta, idname+"_pt_endcaps_"+tname, "pt_PLOT_abseta_bin1_");
                    reflbl = refbk; datalbl = databk;
                } else {
                    extraSpam = "        |#eta| < 1.2"; fits[i][j][0] = single(fit_pt_eta, idname+"_pt_barrel_"+tname,  "pt_PLOT_abseta_bin0_");
                    extraSpam = "  1.2 < |#eta| < 2.4"; fits[i][j][1] = single(fit_pt_eta, idname+"_pt_endcaps_"+tname, "pt_PLOT_abseta_bin1_");
                }
            } else if (ref_pt_eta) {
                refs[i][j][0] = getFit(ref_pt_eta, "pt_PLOT_abseta_bin0_");
                refs[i][j][1] = getFit(ref_pt_eta, "pt_PLOT_abseta_bin1_");
            }

            yMin = 0.0; yMax = 1.1;
            fit_pt_eta = gFile->GetDirectory(basedir+"/"+idname+"_pt_abseta_wide_"+tname+"/");
            ref_pt_eta = ref ? ref->GetDirectory(basedir+"/"+idname+"_pt_abseta_wide_"+tname+"/") : 0;
            if (fit_pt_eta && ref_pt_eta) {
                fits_wide[i][j][0] = getFit(fit_pt_eta, "pt_PLOT_abseta_bin0_");
                fits_wide[i][j][1] = getFit(fit_pt_eta, "pt_PLOT_abseta_bin1_");
                refs_wide[i][j][0] = getFit(ref_pt_eta, "pt_PLOT_abseta_bin0_");
                refs_wide[i][j][1] = getFit(ref_pt_eta, "pt_PLOT_abseta_bin1_");
                extraSpam = "        |#eta| < 1.2"; refstack(fit_pt_eta, ref_pt_eta, idname+"_pt_barrel_wide_"+tname,  "pt_PLOT_abseta_bin0_");
                //extraSpam = "  0.9 < |#eta| < 1.2"; refstack(fit_pt_eta, ref_pt_eta, idname+"_pt_overlap_"+tname, "pt_PLOT_abseta_bin1_");
                extraSpam = "  1.2 < |#eta| < 2.4"; refstack(fit_pt_eta, ref_pt_eta, idname+"_pt_endcaps_wide_"+tname, "pt_PLOT_abseta_bin1_");
            } else if (fit_pt_eta) {
                TDirectory *mc_pt_eta  = gFile->GetDirectory(basedir+"/"+idname+"_pt_abseta_wide_"+tname+"_mcTrue/");
                if (mc_pt_eta) {
                    fits_wide[i][j][0] = getFit(fit_pt_eta, "pt_PLOT_abseta_bin0_");
                    fits_wide[i][j][1] = getFit(fit_pt_eta, "pt_PLOT_abseta_bin1_");
                    TString databk = datalbl, refbk = reflbl;
                    datalbl = "T&P fit"; reflbl = "Sim. truth";
                    extraSpam = "        |#eta| < 1.2"; mcstack(fit_pt_eta, mc_pt_eta, idname+"_pt_barrel_wide_"+tname,  "pt_PLOT_abseta_bin0_");
                    extraSpam = "  1.2 < |#eta| < 2.4"; mcstack(fit_pt_eta, mc_pt_eta, idname+"_pt_endcaps_wide_"+tname, "pt_PLOT_abseta_bin1_");
                    reflbl = refbk; datalbl = databk;
                } else {
                    extraSpam = "        |#eta| < 1.2"; fits_wide[i][j][0] = single(fit_pt_eta, idname+"_pt_barrel_wide_"+tname,  "pt_PLOT_abseta_bin0_");
                    extraSpam = "  1.2 < |#eta| < 2.4"; fits_wide[i][j][1] = single(fit_pt_eta, idname+"_pt_endcaps_wide_"+tname, "pt_PLOT_abseta_bin1_");
                }
            } else if (ref_pt_eta) {
                refs_wide[i][j][0] = getFit(ref_pt_eta, "pt_PLOT_abseta_bin0_");
                refs_wide[i][j][1] = getFit(ref_pt_eta, "pt_PLOT_abseta_bin1_");
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
                        //doCanvas(fit_vtx_barrel,  1, 10, idname+"_vtx_barrel_" +tname+"_vtx_%d",  ".*tag_nVertices_bin%d_");
                        //doCanvas(fit_vtx_endcaps, 1, 10, idname+"_vtx_endcaps_"+tname+"_vtx_%d",  ".*tag_nVertices_bin%d_");
                    }
                }
            }

        }

        yMin = 0.85; yMax = 1.019;
        TDirectory *fit_plateau  = gFile->GetDirectory(basedir+"/"+idname+"_plateau/");
        if (fit_plateau) {
            TDirectory *ref_plateau  = ref ? ref->GetDirectory(basedir+"/"+idname+"_plateau/") : 0;
            if (ref_plateau) {
                extraSpam = (idname == "VBTF" ? "  p_{T} > 10 GeV/c" : "  p_{T} > 6 GeV/c"); 
                refstack(fit_plateau,  ref_plateau,  idname+"_plateau",  "abseta_PLOT_");
            } else {
                TDirectory *mc_plateau  = gFile->GetDirectory(basedir+"/"+idname+"_plateau_mcTrue/");
                if (mc_plateau) {
                    TString databk = datalbl, refbk = reflbl;
                    datalbl = "T&P fit"; reflbl = "Sim. truth";
                    mcstack(fit_plateau,  mc_plateau,  idname+"_plateau",  "abseta_PLOT_");
                    reflbl = refbk; datalbl = databk;
                } else {
                    single(fit_plateau,  idname+"_plateau",  "abseta_PLOT_");
                }
                //doCanvas(fit_plateau,  1, 1, idname+"_plateau" ,  "abseta_bin1_");
            }
        }

    }

    yMin = 0.0; yMax = 1.1;
    bool doFillMCSave = doFillMC;
    doFillMC = false;
    c1->cd(); c1->Clear(); 
    if (doSquare) squareCanvas(c1);
    double xmax = 0, pmax = 0; TH1F *frame = 0, *pframe = 0;
    for (int i = 0; i < nids; ++i) {
        int nfound = 0, jfound = -1;
        for (int j = 0; j < ntrig; ++j) {
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
            for (int j = 0; j < ntrig; ++j) {
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
        for (int i = 0; i < nids; ++i) {
            for (int be = 0; be < neta; ++be) {
                TString label = TString::Format("%s_%s", ids[i], (be == 0 ? "barrel" : (be == 1 ? "endcaps" : "endcaps_p") ));
                TGraphAsymmErrors *mref = merge(&refs[i][0][0], ntrig, neta, be);
                TGraphAsymmErrors *mmct = merge(&mcts[i][0][0], ntrig, neta, be);
                TGraphAsymmErrors *mfit = merge(&fits[i][0][0], ntrig, neta, be);
                if (mfit == 0) { std::cerr << "No mfit for " << label << std::endl; continue; }
                if (mref == 0) { std::cerr << "No mref for " << label << std::endl; continue; }
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
                if (doTxt) printGraphs(mfit, mref, label+"_merge");
                if (doPdf) c1->Print(prefix+label+"_merge.pdf");
                if (fOut) fOut->WriteTObject(mfit, "fit_"+label+"_merge");
                if (fOut) fOut->WriteTObject(mref, "ref_"+label+"_merge");
                if (fOut && mmct) fOut->WriteTObject(mmct, "mct_"+label+"_merge");

                label += "_wide";
                mref = merge(&refs_wide[i][0][0], ntrig, neta, be);
                mmct = merge(&mcts_wide[i][0][0], ntrig, neta, be);
                mfit = merge(&fits_wide[i][0][0], ntrig, neta, be);
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
                if (doTxt) printGraphs(mfit, mref, label+"_merge");
                if (doPdf) c1->Print(prefix+label+"_merge.pdf");
                if (fOut) fOut->WriteTObject(mfit, "fit_"+label+"_merge");
                if (fOut) fOut->WriteTObject(mref, "ref_"+label+"_merge");
                if (fOut && mmct) fOut->WriteTObject(mmct, "mct_"+label+"_merge");

            }
        }
    }
}


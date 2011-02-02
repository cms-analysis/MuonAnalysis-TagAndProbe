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

    doRatioPlot = false;
    doDiffPlot = false;
    doPdf = false;
    doSquare = true; yMin = 0; yMax = 1.1;
    datalbl = "Data, 2010B";
    reflbl  = "Simulation";
    preliminary = ""; //CMS Preliminary,   #sqrt{s} = 7 TeV";
    plotMuonIDData();
}

void plotMuonIDData() {
    retitle = "Efficiency";

    const int nids  = 4;
    char *ids[nids]    = { "Glb",    "TMOST",   "VBTF",  "PF" };
    char *titles[nids] = { "Global",  "Soft",   "Tight", "PF" };

    const int ntrig = 3;
    char *trig[ntrig] = { "Mu5_Track0", "Mu3_Track3", "Mu3_Track5" };
    TGraphAsymmErrors *fits[nids][ntrig][2];
    for (size_t i = 0; i < nids; ++i) { for (size_t j = 0; j < ntrig; ++j) { fits[i][j][0] = 0; fits[i][j][1] = 0; } }
    
    for (size_t i = 0; i < nids; ++i) {
        for (size_t j = 0; j < ntrig; ++j) {
            TString idname(ids[i]);
            TString tname(trig[j]);
            retitle = TString(titles[i])+" muon efficiency";
            TDirectory *fit_pt_eta = gFile->GetDirectory(basedir+"/"+idname+"_pt_abseta_"+tname+"/");
            if (fit_pt_eta == 0) { if (i == 0) { gFile->GetDirectory(basedir)->ls(); } ; continue; }
            if (ref) {
                TDirectory *ref_pt_eta = ref->GetDirectory(basedir+"/"+idname+"_pt_abseta_"+tname+"/");
                if (TString(ref->GetName()) == gFile->GetName()) {
                    ref_pt_eta = ref->GetDirectory(basedir+"/"+idname+"_pt_abseta_sg_"+tname+"/");
                    datalbl = "Data, All"; reflbl = "Data, Seagull";
                }
                extraSpam = "        |#eta| < 1.2"; refstack(fit_pt_eta, ref_pt_eta, idname+"_pt_barrel_"+tname,  "pt_PLOT_abseta_bin0_");
                extraSpam = "  1.2 < |#eta| < 2.4"; refstack(fit_pt_eta, ref_pt_eta, idname+"_pt_endcaps_"+tname, "pt_PLOT_abseta_bin1_");
                TDirectory *mc_pt_eta  = ref->GetDirectory(basedir+"/"+idname+"_pt_abseta_"+tname+"_mcTrue/");
                if (0 && mc_pt_eta) {
                    extraSpam = "             |#eta| < 1.2"; refstack3(fit_pt_eta, ref_pt_eta, mc_pt_eta, idname+"_pt_barrel_3_"+tname,  "pt_PLOT_abseta_bin0_");
                    extraSpam = "       1.2 < |#eta| < 2.4"; refstack3(fit_pt_eta, ref_pt_eta, mc_pt_eta, idname+"_pt_endcaps_3_"+tname, "pt_PLOT_abseta_bin1_");
                }
            } else {
                TDirectory *mc_pt_eta  = gFile->GetDirectory(basedir+"/"+idname+"_pt_abseta_"+tname+"_mcTrue/");
                if (mc_pt_eta) {
                    datalbl = "T&P fit"; reflbl = "Sim. truth";
                    extraSpam = "        |#eta| < 1.2"; mcstack(fit_pt_eta, mc_pt_eta, idname+"_pt_barrel_"+tname,  "pt_PLOT_abseta_bin0_");
                    extraSpam = "  1.2 < |#eta| < 2.4"; mcstack(fit_pt_eta, mc_pt_eta, idname+"_pt_endcaps_"+tname, "pt_PLOT_abseta_bin1_");
                } else {
                    extraSpam = "        |#eta| < 1.2"; fits[i][j][0] = single(fit_pt_eta, idname+"_pt_barrel_"+tname,  "pt_PLOT_abseta_bin0_");
                    extraSpam = "  1.2 < |#eta| < 2.4"; fits[i][j][1] = single(fit_pt_eta, idname+"_pt_endcaps_"+tname, "pt_PLOT_abseta_bin1_");
                }
            }
            if (ref == 0) {
                doCanvas(fit_pt_eta, 1, 37, idname+"_barrel_"+tname+"_pt_%d",   "abseta_bin0__pt_bin%d_");
                doCanvas(fit_pt_eta, 1, 37, idname+"_endcaps_"+tname+"_pt_%d",  "abseta_bin1__pt_bin%d_");
            }
        }
    }

    c1->cd(); c1->Clear(); 
    if (doSquare) squareCanvas(c1);
    for (size_t i = 0; i < nids; ++i) {
        double xmax = 0; int nfound = 0, jfound = -1;
        for (size_t j = 0; j < ntrig; ++j) {
            if (fits[i][j][0] && fits[i][j][1]) {
                xmax = TMath::Max(xmax, xmaxGraph(fits[i][j][0])); 
                xmax = TMath::Max(xmax, xmaxGraph(fits[i][j][1])); 
                jfound = j;
            }
        }
        if (jfound == -1) continue;
        TH1F frame("frame","frame",1,0.,xmax); 
        frame.GetYaxis()->SetRangeUser(yMin,yMax);
        frame.GetYaxis()->SetTitle(fits[i][jfound][0]->GetYaxis()->GetTitle());
        frame.GetXaxis()->SetTitle(fits[i][jfound][0]->GetXaxis()->GetTitle());
        int founds[ntrig];
        int colors[6] = { 1, 2, 4, 6, 9, 206 };
        for (int be = 0; be <= 1; ++be) {
            frame.Draw();
            extraSpam = (be == 0 ? "        |#eta| < 1.2" : "  1.2 < |#eta| < 2.4");
            nfound = 0;
            for (size_t j = 0; j < ntrig; ++j) {
                if (fits[i][j][be]) { 
                    fits[i][j][be]->SetLineColor(colors[nfound]);
                    fits[i][j][be]->SetMarkerColor(colors[nfound]);
                    if (ref == 0) fits[i][j][be]->Draw("P");
                    founds[nfound++] = j;
                }
            }
            if (ref == 0) {
                if (nfound == 2) doLegend(fits[i][founds[0]][be], fits[i][founds[1]][be], trig[founds[0]], trig[founds[1]]);
                else if (nfound == 3) doLegend(fits[i][founds[0]][be], fits[i][founds[1]][be], fits[i][founds[2]][be], trig[founds[0]], trig[founds[1]], trig[founds[1]]);
                TString label = TString::Format("%s_%s", ids[i], (be == 0 ? "barrel" : "endcaps"));
                c1->Print(prefix+label+"_all.png");
                if (doPdf) c1->Print(prefix+label+"_all.pdf");
            }
            std::cout << "found: " << nfound << std::endl;
            if (nfound > 1) {
                TGraphAsymmErrors *g1 = fits[i][founds[0]][be], *g2 = fits[i][founds[1]][be];
                TGraphAsymmErrors *g3 = nfound > 2 ? fits[i][founds[2]][be] : NULL;
                TGraphAsymmErrors *g4 = nfound > 3 ? fits[i][founds[3]][be] : NULL;
                TGraphAsymmErrors *merged = mergeGraphs(g1,g2,g3,g4);
                frame.Draw();
                merged->Draw("P");
                c1->Print(prefix+label+"_merge.png");
                if (doPdf) c1->Print(prefix+label+"_merge.pdf");
            }
        }
    }
}


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
TString prefix = "plots_dev/tracking/";
TString basedir  = "tpTreeSta";

TFile *ref = 0;

TCanvas *c1 = 0;
void plotTracking(TString scenario="data",TString match="dr050e020") {
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
    doDiffPlot  = false;
    doPdf = false;
    doSquare = true; yMin = 0.8; yMax = 1.1;
    datalbl = "Data";
    reflbl  = "Sim.";
    preliminary = ""; //CMS Preliminary, #sqrt{s} = 7 TeV";
    plotTracking_(match);
}

void plotTracking_(TString match) {

    const int nplots = 2;
    const char *plots[nplots] = { "eff_ichep", "eff_vxt" };
    const char * vars[nplots] = { "abseta",    "tag_nVertices" };
    for (size_t i = 0; i < nplots; ++i) {
        TString plotname(plots[i]); plotname += "_"+match;
        TString varname(vars[i]);
        TDirectory *fit_0 = gFile->GetDirectory(basedir+"/"+plotname+"/");
        if (fit_0 == 0) { std::cerr << "Didn't find " << basedir+"/"+plotname+"/" << " in " << gFile->GetName() << std::endl; continue; }
        TDirectory *fit_1 = gFile->GetDirectory(basedir+"/"+plotname+"NoJPsi/");
        TDirectory *fit_2 = gFile->GetDirectory(basedir+"/"+plotname+"NoBestJPsi/");
        if (ref != 0) {
            TDirectory *ref_0 = ref->GetDirectory(basedir+"/"+plotname+"/");
            TDirectory *ref_1 = ref->GetDirectory(basedir+"/"+plotname+"NoJPsi/");
            TDirectory *ref_2 = ref->GetDirectory(basedir+"/"+plotname+"NoBestJPsi/");
            yMin = 0.9; yMax = 1.05;
            retitle = "Raw efficiency";
            if (fit_0 && ref_0) refstack(fit_0, ref_0, plotname, varname+"_PLOT_");
            retitle = "Fake rate";
            yMin = 0.0; yMax = 1.1;
            if (fit_1 && ref_1) refstack(fit_1, ref_1, plotname+"_fake1", varname+"_PLOT_");
            if (fit_2 && ref_2) refstack(fit_2, ref_2, plotname+"_fake2", varname+"_PLOT_");
        } else {
            TDirectory *mc_pt_eta = 0;
            if (mc_pt_eta) {
                //datalbl = "T&P fit"; reflbl = "Sim. truth";
                //extraSpam = "        |#eta| < 1.2"; mcstack(fit_pt_eta, mc_pt_eta, idname+"_pt_barrel",  "pt_PLOT_abseta_bin0_");
                //extraSpam = "  1.2 < |#eta| < 2.4"; mcstack(fit_pt_eta, mc_pt_eta, idname+"_pt_endcaps", "pt_PLOT_abseta_bin1_");
            } else {
                yMin = 0.9; yMax = 1.05; retitle = "Raw efficiency";
                if (fit_0) single(fit_0, plotname, varname+"_PLOT_");
                yMin = 0.0; yMax = 1.1;  retitle = "Fake rate";
                if (fit_1) single(fit_1, plotname+"_fake1", varname+"_PLOT_");
                if (fit_2) single(fit_2, plotname+"_fake2", varname+"_PLOT_");
                yMin = 0.9; yMax = 1.05; retitle = "Corrected efficiency";
                if (fit_0 && fit_1) corrsingle(fit_0, fit_1, plotname+"_corr1", varname+"_PLOT_");
                if (fit_0 && fit_2) corrsingle(fit_0, fit_2, plotname+"_corr2", varname+"_PLOT_");
            }
        }

        if (1) {
            if (fit_0) doCanvas(fit_0, 1, 50, plotname+"_"+varname+"_%d", varname+"_bin%d__");
            if (fit_1) doCanvas(fit_1, 1, 50, plotname+"_fake1_"+varname+"_%d", varname+"_bin%d__");
            if (fit_2) doCanvas(fit_1, 1, 50, plotname+"_fake2_"+varname+"_%d", varname+"_bin%d__");
        }
    }
}


void corrsingle(TDirectory *fit, TDirectory *fake, TString alias, TString fitname) {
    if (fake == 0 || fit == 0) return;
    TCanvas *pfake = getFromPrefix(fake->GetDirectory("fit_eff_plots"), fitname);
    if (pfake == 0) {
        std::cerr << "NOT FOUND: " << "fit_eff_plots/"+fitname << " in " << fake->GetName() << std::endl;
        return;
    }
    RooHist *hfake = (RooHist *) pfake->FindObject("hxy_fit_eff");

    TCanvas *pfit = getFromPrefix(fit->GetDirectory("fit_eff_plots"), fitname);
    if (pfit == 0) {
        std::cerr << "NOT FOUND: " << "fit_eff_plots/"+fitname << " in " << fit->GetName() << std::endl;
        return;
    }
    RooHist *hfit = (RooHist *) pfit->FindObject("hxy_fit_eff");

    TGraphAsymmErrors *out = new TGraphAsymmErrors();
    for (int i = 0, n = hfit->GetN(), k = 0; i < n; ++i) {
        double x = hfit->GetX()[i], y = hfit->GetY()[i], eyh = hfit->GetErrorYhigh(i), eyl = hfit->GetErrorYlow(i);
        int j = findBin(hfake, x); if (j == -1) continue;
        double yf  = hfake->GetY()[j], eyhf = hfake->GetErrorYhigh(j), eylf = hfake->GetErrorYlow(j);
        double ycorr =    (   y   -    yf   )/(1-    yf   );
        double ycorr_hi = ((y+eyh)-(yf-eylf))/(1-(yf-eylf));
        double ycorr_lo = ((y-eyl)-(yf+eyhf))/(1-(yf+eyhf));
        /*
        std::cout << "x = " << x << " [" << (x-hfit->GetErrorXlow(i)) << ", " << (x+hfit->GetErrorXhigh(i)) << "] \t" <<
                     "y  = " << y  << " -"<<eyl <<"/+"<<eyh << " \t " <<
                     "yf = " << yf << " -"<<eylf<<"/+"<<eyhf<< " \t " <<
                     "yc = " << ycorr << "[ " << ycorr_lo << ", " << ycorr_hi << "]" << std::endl;
        */
        out->Set(k+1);
        out->SetPoint(k, x, ycorr);
        out->SetPointError(k, hfit->GetErrorXlow(i), hfit->GetErrorXhigh(i), ycorr - ycorr_lo, ycorr_hi - ycorr );
        k++;
    }

    c1->Clear(); c1->cd();
    out->SetLineWidth(2);
    out->SetLineColor(kBlack);
    out->SetMarkerColor(kBlack);
    out->SetMarkerStyle(20);
    out->SetMarkerSize(1.6);
    out->Draw("AP");
    out->GetXaxis()->SetTitle(getXtitle(pfit)); 
    out->GetXaxis()->SetMoreLogLabels(1);
    out->GetXaxis()->SetRangeUser(out->GetX()[0]-out->GetErrorXlow(0), out->GetX()[out->GetN()-1]+out->GetErrorXhigh(out->GetN()-1));
    out->GetYaxis()->SetRangeUser(yMin, yMax);
    if (retitle != "") out->GetYaxis()->SetTitle(retitle);
    
    if (doSquare) squareCanvas(c1);
    if (preliminary != "") cmsprelim();

    gPad->Print(prefix+alias+".png");
    if (doPdf) gPad->Print(prefix+alias+".pdf");
}


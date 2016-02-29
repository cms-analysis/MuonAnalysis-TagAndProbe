/**
  root.exe -b -l -q data.root data_NoJPsi.root mc.root mc_NoJPsi.root
*/
#include <TCanvas.h>
#include <TPad.h>
#include "plotUtil_Tracking.cxx"

TString basedir  = "tpTree";
TString basedir2 = "tpTree";

TFile *files[10]; int gNFiles = 0;
void plotTracking_(TString match) ;
TGraphAsymmErrors* corrsingle(TDirectory *fit, TDirectory *fake, TString alias, TString fitname, bool print=true) ;
TGraphAsymmErrors* corrsingle2(TDirectory *fit, TDirectory *fake, TDirectory *fake2, TString alias, TString fitname, bool print=true) ;
TGraphAsymmErrors* biasCorrection(TGraphAsymmErrors *fit0, TGraphAsymmErrors *fit1, TGraphAsymmErrors *ratioSP) ;

void plotTracking(TString scenario="data",bool DataVsMCcomparison = false, TString match="dr030e030") {
    prefix = prefix+scenario+"/";
    basedir  = "tpTreeSta";
    basedir2 = "tpTreeSta";

    gROOT->ProcessLine(".x /afs/cern.ch/user/g/gpetrucc/cpp/tdrstyle.cc");
    gStyle->SetOptStat(0);
    c1 = new TCanvas("c1","c1");

    if (gROOT->GetListOfFiles()->GetEntries() >= 4) {
        for (int i = 0, n = gROOT->GetListOfFiles()->GetEntries(); i < n; ++i) {
            files[i] = (TFile*) gROOT->GetListOfFiles()->At(i);
        }
        gNFiles = gROOT->GetListOfFiles()->GetEntries();
    } else return;

    doRatioPlot = false;
    doDiffPlot  = false;
    doPdf = true;
    doSquare = true; yMin = 0.949; yMax = 1.009;
    doFillMC = true;
    datalbl = "Data, 2015";
    reflbl  = "Simulation";
//    cmslabel = "CMS Simulation, #sqrt{s} = 13 TeV";

    if (scenario.Contains("74X")) {
        datalbl = "74X Data";
        reflbl  = "53X Data";
        lineColorRef_FillMC = 100;
        fillColorRef_FillMC = 207;
        cmslabel = "Run2012D";
        prefix = "plots/";
    }
    if (scenario.Contains("_vs_")) {
        datalbl = "Modified";
        reflbl  = "Reference";
        prefix = "plots/"+scenario+"/";
        yMin = 0.9;
        if (scenario.Contains("_L1_vs_")) {
            basedir  = "tpTreeL1";
            basedir2 = "tpTreeL1";
        }
        if (scenario.Contains("_vs_MCGen")) {
            datalbl = "Tag & Probe";
            reflbl  = "Tag & Gen ";
            basedir2 = "tpTreeGen";
            lineColorRef_FillMC = kOrange+10;
            fillColorRef_FillMC = kOrange-3;
            if (scenario.Contains("_L1_vs_")) {
                datalbl = "Tag & L1";
            }
        }
        if (DataVsMCcomparison == true) {
            cmslabel = "CMS Preliminary";
  	    cmsconditionslabel = "#sqrt{s} = 13 TeV, L = 2.23 fb^{-1}";
            datalbl = "Data, 2015";
            reflbl = "Simulation";
            yMin = 0.9;
        }

    }

    gSystem->mkdir(prefix,true);
    gSystem->Exec("mkdir -p "+prefix);
    gSystem->Exec("cp /afs/cern.ch/user/g/gpetrucc/php/index.php "+prefix+"/");
    fOut = TFile::Open(prefix+"/fits.root","UPDATE");
    plotTracking_(match);
    prefix += "tk0/";
    gSystem->Exec("mkdir -p "+prefix);
    gSystem->Exec("cp /afs/cern.ch/user/g/gpetrucc/php/index.php "+prefix+"/");
    fOut = TFile::Open(prefix+"/fits.root","UPDATE");
    if (scenario.Contains("_vs_")) {
        lineColorRef_FillMC =  kViolet-2;
        fillColorRef_FillMC =  kViolet-2;
    }
    plotTracking_("tk0_"+match);
}

void plotTracking_(TString match) {

  const int nplots = 14;
  const char *plots[nplots] = { "eff_abseta",    "eff_1",    "efft_1",    "eff_eta",  "eff_eta2", "eff_eta3",  "eff",       "eff_vtx",             "eff_two", 	"eff_phi", 	"eff_eph", 	"eff_eph", 	"eff_eph", "eff_pt"};
  const char * vars[nplots] = { "abseta",      "eta",      "eta",       "eta",      "eta",      "eta",       "eta",       "tag_nVertices",       "abseta", 	"phi",	 "phi", 	"phi", 		"phi", "pt"};
  const char *xvars[nplots] = { "muon |#eta|", "muon #eta","muon #eta", "muon #eta","muon #eta","muon #eta", "muon #eta", "N(primary vertices)", "muon |#eta|", "muon #phi",	"muon #phi", 	"muon #phi", 	"muon #phi", "pt"};
  const char *bincs[nplots] = { "",            "",         "",          "",         "",         "",          "",          "",                    "" , 	"",		"eta_bin0", 	"eta_bin1", 	"eta_bin2", ""}; 
  const char *binls[nplots] = { "",            "",         "",          "",         "",         "",          "",          "",                    "" , 	"", 		"#eta < -0.5", 	"-0.5 > #eta > 0.5", "#eta > 0.5", ""};
  bool  has10files = (gNFiles == 10);
  bool  has9files = (gNFiles == 9);
  bool  has8files = (gNFiles >= 8);
  bool  has6files = (gNFiles == 6);
  bool  has4files = (gNFiles == 4) || has8files;
std::cout << "start plotTracking with " << gNFiles << " files." << std::endl;
  TFile *ref = (TFile *) gROOT->GetListOfFiles()->At(1); 
  for (size_t i = 0; i < nplots; ++i) {
    retitleX = TString(xvars[i]);
    TString plotName(plots[i]); plotName += "_"+match;
    TString varname(vars[i]);
    TString binCut(bincs[i]); 
    TString binName = (binCut != "" ? binCut : binCut);
    extraSpam = TString(binls[i]);
    TDirectory *fit_0, *fit_1, *fit_2;
	
    if (has4files || has6files) files[0]->cd();
    fit_0 = gFile->GetDirectory(basedir+"/"+plotName+"/");
    if (fit_0 == 0) { 
      std::cerr << "Didn't find " << basedir+"/"+plotName+"/" << " in " << gFile->GetName() << std::endl; 
      continue; 
    }

    if (has4files || has6files) files[1]->cd();
    fit_1 = gFile->GetDirectory(basedir+"/"+plotName+"NoJPsi/");
    if (fit_1 == 0) {
      std::cerr << "Didn't find " << basedir+"/"+plotName+"/" << " in " << gFile->GetName() << std::endl;
      continue;
    }

    if (has6files) {
      files[2]->cd();
      fit_2 = gFile->GetDirectory(basedir+"/"+plotName+"NoBestJPsi/");
      if (fit_2 == 0) {
        std::cerr << "Didn't find " << basedir+"/"+plotName+"/" << " in " << gFile->GetName() << std::endl;
        continue;
      }
    }

    if (ref != 0 || has4files || has6files || has8files) {
std::cout << "here.." << std::endl;

            //definition of reference files
            TDirectory *ref_0, *ref_1, *ref_2;
            if (has4files) files[2]->cd();//ref = files[2];
            if (has6files) files[3]->cd();//ref = files[3];
            ref_0 = gFile->GetDirectory(basedir2+"/"+plotName+"/");
            if (has4files) files[3]->cd();//ref = files[3];
            if (has6files) files[4]->cd();//ref = files[4];
            ref_1 = gFile->GetDirectory(basedir2+"/"+plotName+"NoJPsi/");
            if (has6files) files[5]->cd();//ref = files[4];
            ref_2 = gFile->GetDirectory(basedir2+"/"+plotName+"NoBestJPsi/");

            //refstack for eff and fake
            double yMin0 = yMin;; yMax = 1.009;
            retitle = "Raw efficiency";
            if (fit_0 && ref_0) refstack(fit_0, ref_0, plotName+binName, varname+"_PLOT_"+binCut);
            retitle = "Fake match rate";
            yMin = 0.0; yMax = 1.1;
            if (fit_1 && ref_1) refstack(fit_1, ref_1, plotName+binName+"_fake1", varname+"_PLOT_"+binCut);
            if (fit_2 && ref_2) refstack(fit_2, ref_2, plotName+binName+"_fake2", varname+"_PLOT_"+binCut);
            yMin = yMin0;

            if (fit_0 && fit_1 && ref_0 && ref_1) {
std::cout << "Filling graphs 0&1.." << std::endl;
                yMin = yMin0; yMax = 1.003; retitle = "Efficiency";
                TGraphAsymmErrors *corr, *cref;
                if(has4files){
   		  corr = corrsingle(fit_0, fit_1, plotName+"_corr",     varname+"_PLOT_"+binCut, false);
                  cref = corrsingle(ref_0, ref_1, plotName+"_corr_ref", varname+"_PLOT_"+binCut, false);
                } else if ( fit_2 && ref_2 && has6files ){
std::cout << "Filling graphs 2.." << std::endl;
   		  corr = corrsingle2(fit_0, fit_1, fit_2, plotName+"_corr",     varname+"_PLOT_"+binCut, false);
                  cref = corrsingle2(ref_0, ref_1, ref_2, plotName+"_corr_ref", varname+"_PLOT_"+binCut, false);
		}
                if (doFillMC) {
                    cref->SetLineColor(lineColorRef_FillMC);
                    cref->SetFillColor(fillColorRef_FillMC);
                    cref->SetLineStyle(0);
                    cref->SetMarkerColor(2);
                    cref->SetMarkerStyle(21);
                    cref->SetMarkerSize(0.4);
                } else {
                    cref->SetLineWidth(2);
                    cref->SetLineColor(kRed);
                    cref->SetMarkerColor(kRed);
                    cref->SetMarkerStyle(25);
                    cref->SetMarkerSize(2.0);
                }
                corr->SetLineWidth(2);
                corr->SetLineColor(kBlack);
                corr->SetMarkerColor(kBlack);
                corr->SetMarkerStyle(20);
                corr->SetMarkerSize(1.6);
                cref->GetXaxis()->SetTitle(retitleX);
                cref->GetYaxis()->SetTitleOffset(1.6);
                double myYMin = TMath::Min(yminGraph(corr),yminGraph(cref));
                myYMin = TMath::Max(0.0, 1.0 - 1.45*(1.0-myYMin));
                cref->GetYaxis()->SetRangeUser(yMin > yMax ? myYMin : yMin, yMax);
                gStyle->SetErrorX(0.5);
                gPad->SetLeftMargin(0.2);
                cref->Draw(doFillMC ? "AE2" : "AP");
                corr->Draw("P SAME");
                if (datalbl) doLegend(corr,cref,datalbl,reflbl);
                c1->Print(prefix+plotName+binName+"_corr"+".png");
                if (fOut) { fOut->WriteTObject(corr,plotName+binName+"_corr", "Overwrite"); }
                if (fOut) { fOut->WriteTObject(cref,"ref_"+plotName+binName+"_corr", "Overwrite"); }
                if (doPdf) c1->Print(prefix+plotName+binName+"_corr"+".pdf");
                if (doTxt)  printGraph(corr,"fit_"+plotName+binName+"_corr");
                if (doTxt)  printGraph(cref,"ref_"+plotName+binName+"_corr");
                if (doTxt)  printGraphs(corr, cref, plotName+binName+"_corr");
                //autoScale = (yMin >= 0); yMinR = 0.964; yMaxR = 1.007;
                TString bk_datalbl = datalbl, bk_reflbl = reflbl;
                datalbl = "Data"; reflbl = "Sim.";
                if (bk_datalbl == "Modified") datalbl = "Mod.";
                if (bk_reflbl == "Reference") reflbl = "Ref.";
                doRatio(corr,cref,plotName+binName+"_corr",retitleX); 
                reflbl = bk_reflbl; datalbl = bk_datalbl;

                if (has8files) {
                    TString plotName0(plots[i]); plotName0 += "_tk0_"+match;
                    TDirectory *fit0_0 = files[4]->GetDirectory(basedir+"/"+plotName0+"/");
                    TDirectory *fit0_1 = files[5]->GetDirectory(basedir+"/"+plotName0+"NoJPsi/");
                    TDirectory *ref0_0 = files[6]->GetDirectory(basedir2+"/"+plotName0+"/");
                    TDirectory *ref0_1 = files[7]->GetDirectory(basedir2+"/"+plotName0+"NoJPsi/");
                    if (!fit0_0) std::cout << "Missing fit0_0" << std::endl;
                    if (!fit0_1) std::cout << "Missing fit0_1" << std::endl;
                    if (!ref0_0) std::cout << "Missing ref0_0" << std::endl;
                    if (!ref0_1) std::cout << "Missing ref0_1" << std::endl;
                    TGraphAsymmErrors *corr0 = corrsingle(fit0_0, fit0_1, plotName0+"_corr0",     varname+"_PLOT_"+binCut, false);
                    TGraphAsymmErrors *cref0 = corrsingle(ref0_0, ref0_1, plotName0+"_corr0_ref", varname+"_PLOT_"+binCut, false);
                    cref0->SetFillColor(kViolet-2);
                    cref0->SetLineStyle(0);
                    corr0->SetFillColor(kViolet+3);
                    corr0->SetLineColor(kViolet+3);
                    corr0->SetLineWidth(2);
                    corr0->SetMarkerStyle(24);
                    corr0->SetMarkerSize(1.3);
                    myYMin = TMath::Min(yminGraph(corr),yminGraph(cref));
                    myYMin = TMath::Min(myYMin, TMath::Min(yminGraph(corr0),yminGraph(cref0)));
                    myYMin = TMath::Max(0.0, 1.0 - 1.45*(1.0-myYMin));
                    cref0->GetYaxis()->SetRangeUser(yMin > yMax ? myYMin : yMin, yMax);
                    cref0->GetXaxis()->SetTitle(retitleX);
                    cref0->Draw(doFillMC ? "AE2" : "AP");
                    cref->Draw(doFillMC ? "E2" : "P");
                    corr0->Draw("P SAME");
                    corr->Draw("P SAME");
                    if (datalbl) doLegend(corr,cref,datalbl,reflbl);
                    c1->Print(prefix+plotName+binName+"_corr4"+".png"); 
                    if (doPdf) c1->Print(prefix+plotName+binName+"_corr4"+".pdf");

                    if (has9files && fOut != 0) {
                        TString plotNameS( plots[i]+3); plotNameS =   "effS"+plotNameS;
                        TString plotNameSP(plots[i]+3); plotNameSP = "effSP"+plotNameSP;
                        plotNameS.ReplaceAll("t_1",""); plotNameSP.ReplaceAll("t_1","");
                        plotNameS.ReplaceAll("_1",""); plotNameSP.ReplaceAll("_1","");
                        TString bk_datalbl = datalbl, bk_reflbl = reflbl;
                        retitle = "OITK Seed Efficiency";
                        datalbl = "Tk Probe"; reflbl = "Tk & Sta";
                        TDirectory *fitS  = files[8]->GetDirectory("tpTree/"+plotNameS+"/");
                        TDirectory *fitSP = files[8]->GetDirectory("tpTree/"+plotNameSP+"/");
                        doRatioPlot = 1;
                        if (fitS && fitSP) refstack(fitS, fitSP, plotNameS+binName, varname+"_PLOT_"+binCut);
                        reflbl = bk_reflbl; datalbl = bk_datalbl;
                        TString corrName = "ref_"+plotNameS+binName;
                        if (*(plots[i]+3) == 't') corrName = "fit_"+plotNameS+binName;
                        std::cout << "Correction using " << corrName << std::endl;
                        TGraphAsymmErrors *ratioSP = (TGraphAsymmErrors *) fOut->Get(corrName);
                        TGraphAsymmErrors *corrBC = biasCorrection(corr0, corr, ratioSP);
                        corrBC->SetLineWidth(2);
                        corrBC->SetLineColor(kBlack);
                        corrBC->SetMarkerColor(kBlack);
                        corrBC->SetMarkerStyle(20);
                        corrBC->SetMarkerSize(1.6);
                        myYMin = TMath::Min(yminGraph(corrBC),yminGraph(cref));
                        myYMin = TMath::Max(0.0, 1.0 - 1.45*(1.0-myYMin));
                        cref->GetYaxis()->SetRangeUser(yMin > yMax ? myYMin : yMin, yMax);
                        cref->GetYaxis()->SetTitle("Bias-corrected eff.");
                        cref->Draw(doFillMC ? "AE2" : "AP");
                        corrBC->Draw("P SAME");
                        if (datalbl) doLegend(corrBC,cref,datalbl,reflbl);
                        c1->Print(prefix+plotName+binName+"_biascorr"+".png"); 
                        if (doPdf) c1->Print(prefix+plotName+binName+"_biascorr"+".pdf");
                        if (doTxt)  printGraph(corrBC,"fit_"+plotName+binName+"_biascorr");
                        if (doTxt)  printGraphs(corrBC, cref, plotName+binName+"_biascorr");
                        doRatio(corrBC,cref,plotName+binName+"_biascorr",retitleX); 
                    } else if (has10files && fOut != 0) {
                        TString plotNameS( plots[i]+3); plotNameS =   "effS"+plotNameS;
                        TString plotNameSP(plots[i]+3); plotNameSP = "effSP"+plotNameSP;
                        TDirectory *fitS  = files[8]->GetDirectory("tpTree/"+plotNameS+"/");
                        TDirectory *fitSP = files[8]->GetDirectory("tpTree/"+plotNameSP+"/");
                        TDirectory *refS  = files[9]->GetDirectory("tpTree/"+plotNameS+"/");
                        TDirectory *refSP = files[9]->GetDirectory("tpTree/"+plotNameSP+"/");
                        retitle = "OITK Seed Efficiency";
                        if (fitS  && refS)  refstack(fitS,  refS,  plotNameS +binName, varname+"_PLOT_"+binCut);
                        retitle = "OITK Probe Efficiency";
                        if (fitSP && refSP) refstack(fitSP, refSP, plotNameSP+binName, varname+"_PLOT_"+binCut);
                        TString bk_datalbl = datalbl, bk_reflbl = reflbl;
                        datalbl = "Tk Probe"; reflbl = "Tk & Sta";
                        doRatioPlot = 1;
                        if (fitS && fitSP) refstack(fitS, fitSP, plotNameS+binName+"_fit", varname+"_PLOT_"+binCut);
                        if (refS && refSP) refstack(refS, refSP, plotNameS+binName+"_ref", varname+"_PLOT_"+binCut);
                        reflbl = bk_reflbl; datalbl = bk_datalbl;
                        TGraphAsymmErrors *fitRSP = (TGraphAsymmErrors *) fOut->Get("ratio_"+plotNameS+binName+"_fit");
                        TGraphAsymmErrors *refRSP = (TGraphAsymmErrors *) fOut->Get("ratio_"+plotNameS+binName+"_ref");
                        TGraphAsymmErrors *corrBC = biasCorrection(corr0, corr, fitRSP);
                        TGraphAsymmErrors *crefBC = biasCorrection(cref0, cref, refRSP);
                        corrBC->SetLineWidth(2);
                        corrBC->SetLineColor(kBlack);
                        corrBC->SetMarkerColor(kBlack);
                        corrBC->SetMarkerStyle(20);
                        corrBC->SetMarkerSize(1.6);
                        crefBC->SetLineColor(lineColorRef_FillMC);
                        crefBC->SetFillColor(fillColorRef_FillMC);
                        crefBC->SetLineStyle(0);
                        crefBC->SetMarkerColor(2);
                        crefBC->SetMarkerStyle(21);
                        crefBC->SetMarkerSize(0.4);
                        myYMin = TMath::Min(yminGraph(corrBC),yminGraph(crefBC));
                        myYMin = TMath::Max(0.0, 1.0 - 1.45*(1.0-myYMin));
                        crefBC->GetYaxis()->SetRangeUser(yMin > yMax ? myYMin : yMin, yMax);
                        crefBC->GetYaxis()->SetTitle("Bias-corrected eff.");
                        crefBC->Draw(doFillMC ? "AE2" : "AP");
                        corrBC->Draw("P SAME");
                        if (datalbl) doLegend(corrBC,crefBC,datalbl,reflbl);
                        c1->Print(prefix+plotName+binName+"_biascorr"+".png"); 
                        if (doPdf) c1->Print(prefix+plotName+binName+"_biascorr"+".pdf");
                        if (doTxt)  printGraph(corrBC,"fit_"+plotName+binName+"_biascorr");
                        if (doTxt)  printGraphs(corrBC, crefBC, plotName+binName+"_biascorr");
                        doRatio(corrBC,cref,plotName+binName+"_biascorr",retitleX); 
                    }
		}                
            }
 
/*
            else {
              TDirectory *mc_pt_eta = 0;
              if (mc_pt_eta) {
                //datalbl = "T&P fit"; reflbl = "Sim. truth";
                //extraSpam = "        |#eta| < 1.2"; mcstack(fit_pt_eta, mc_pt_eta, idname+"_pt_barrel",  "pt_PLOT_abseta_bin0_");
                //extraSpam = "  1.2 < |#eta| < 2.4"; mcstack(fit_pt_eta, mc_pt_eta, idname+"_pt_endcaps", "pt_PLOT_abseta_bin1_");
              } else {
                yMin = 0.9; yMax = 1.019; retitle = "Raw efficiency";
                if (fit_0) single(fit_0, plotName, varname+"_PLOT_"+binCut);
                yMin = 0.0; yMax = 1.1;  retitle = "Fake rate";
                if (fit_1) single(fit_1, plotName+"_fake1", varname+"_PLOT_"+binCut);
                if (fit_2) single(fit_2, plotName+"_fake2", varname+"_PLOT_"+binCut);
                yMin = 0.9; yMax = 1.019; retitle = "Efficiency";
                if (fit_0 && fit_1) corrsingle(fit_0, fit_1, plotName+"_corr1", varname+"_PLOT_"+binCut);
	        if (fit_0 && fit_2) corrsingle(fit_0, fit_2, plotname+"_corr2", varname+"_PLOT_");
		if (fit_0 && fit_1 && fit_2) corrsingle2(fit_0, fit_1, fit_2, plotname+"_corr", varname+"_PLOT_");
              }
          }
        }
      if (0 && ref == 0) {
        if (fit_0) doCanvas(fit_0, 1, 50, plotName+"_"+varname+"_%d", varname+"_bin%d__");
        if (fit_1) doCanvas(fit_1, 1, 50, plotName+"_fake1_"+varname+"_%d", varname+"_bin%d__");
        if (fit_2) doCanvas(fit_2, 1, 50, plotname+"_fake2_"+varname+"_%d", varname+"_bin%d__");
      }
*/
    }
  }//loop over plots

}


TGraphAsymmErrors* corrsingle(TDirectory *fit, TDirectory *fake, TString alias, TString fitname, bool print=true) {
    if (fake == 0 || fit == 0) return 0;
    TCanvas *pfake = (TCanvas*) getFromPrefix(fake->GetDirectory("fit_eff_plots"), fitname);
    if (pfake == 0) {
        std::cerr << "NOT FOUND: " << "fit_eff_plots/"+fitname << " in " << fake->GetName() << std::endl;
        return 0;
    }
    RooHist *hfake = (RooHist *) pfake->FindObject("hxy_fit_eff");

    TCanvas *pfit = (TCanvas*) getFromPrefix(fit->GetDirectory("fit_eff_plots"), fitname);
    if (pfit == 0) {
        std::cerr << "NOT FOUND: " << "fit_eff_plots/"+fitname << " in " << fit->GetName() << std::endl;
        return 0;
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
    out->GetYaxis()->SetDecimals(true);
    out->GetYaxis()->SetTitleOffset(1.3);
    if (retitle != "") out->GetYaxis()->SetTitle(retitle);
    
    if (doSquare) squareCanvas(c1);
    if (cmslabel != "") printcmsprelim();
    if (cmsconditionslabel != "") printcmsconditionslabel();

    if (print) {
        c1->Print(prefix+alias+".png");
        if (doPdf) c1->Print(prefix+alias+".pdf");
    }

    return out;
}

TGraphAsymmErrors* corrsingle2(TDirectory *fit, TDirectory *fake, TDirectory *fake2, TString alias, TString fitname, bool print=true) {
    if (fake == 0 || fake2 == 0 || fit == 0) return 0;
    
    TCanvas *pfake = (TCanvas*) getFromPrefix(fake->GetDirectory("fit_eff_plots"), fitname);
    if (pfake == 0) {
        std::cerr << "NOT FOUND: " << "fit_eff_plots/"+fitname << " in " << fake->GetName() << std::endl;
        return 0;
    }
    RooHist *hfake = (RooHist *) pfake->FindObject("hxy_fit_eff");

    TCanvas *pfake2 = (TCanvas*) getFromPrefix(fake2->GetDirectory("fit_eff_plots"), fitname);
    if (pfake2 == 0) {
        std::cerr << "NOT FOUND: " << "fit_eff_plots/"+fitname << " in " << fake2->GetName() << std::endl;
        return 0;
    }
    RooHist *hfake2 = (RooHist *) pfake2->FindObject("hxy_fit_eff");

    TCanvas *pfit = (TCanvas*) getFromPrefix(fit->GetDirectory("fit_eff_plots"), fitname);
    if (pfit == 0) {
        std::cerr << "NOT FOUND: " << "fit_eff_plots/"+fitname << " in " << fit->GetName() << std::endl;
        return 0;
    }
    RooHist *hfit = (RooHist *) pfit->FindObject("hxy_fit_eff");

    TGraphAsymmErrors *out = new TGraphAsymmErrors();
    for (int i = 0, n = hfit->GetN(), k = 0; i < n; ++i) {
        double x = hfit->GetX()[i], y = hfit->GetY()[i], eyh = hfit->GetErrorYhigh(i), eyl = hfit->GetErrorYlow(i);
        int j = findBin(hfake, x); if (j == -1) continue;
        double yf1  = hfake->GetY()[j], eyhf1 = hfake->GetErrorYhigh(j), eylf1 = hfake->GetErrorYlow(j);
        int j2 = findBin(hfake2, x); if (j2 == -1) continue;
        double yf2  = hfake2->GetY()[j2], eyhf2 = hfake2->GetErrorYhigh(j2), eylf2 = hfake2->GetErrorYlow(j2);
        double yf = 0.5*(yf1+yf2);
        double ylf = TMath::Min(yf1-eylf1, yf2-eylf2);
        double yhf = TMath::Max(yf1+eylf1, yf2+eylf2);
        double ycorr =    (   y   -yf )/(1-yf );
        double ycorr_hi = ((y+eyh)-ylf)/(1-ylf);
        double ycorr_lo = ((y-eyl)-yhf)/(1-yhf);

        out->Set(k+1);
        out->SetPoint(k, x, ycorr);
        out->SetPointError(k, hfit->GetErrorXlow(i), hfit->GetErrorXhigh(i), ycorr - ycorr_lo, ycorr_hi - ycorr );
        k++;
    }

    c1->Clear(); c1->cd();
    c1->SetLeftMargin(0.16);
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
    out->GetYaxis()->SetDecimals(true);
    out->GetYaxis()->SetTitleOffset(1.3);
    if (retitle != "") out->GetYaxis()->SetTitle(retitle);

    if (doSquare) squareCanvas(c1);
    if (cmslabel != "") printcmsprelim();
    if (cmsconditionslabel != "") printcmsconditionslabel();

    if (print) {
        c1->Print(prefix+alias+".png");
        if (doPdf) c1->Print(prefix+alias+".pdf");
    }

    return out;
}

TGraphAsymmErrors* biasCorrection(TGraphAsymmErrors *fit0, TGraphAsymmErrors *fit1, TGraphAsymmErrors *corrFactor) {
    TGraphAsymmErrors *out = new TGraphAsymmErrors();
    for (int i = 0, n = fit0->GetN(), k = 0; i < n; ++i) {
        double x = fit0->GetX()[i], y0 = fit0->GetY()[i], eyh0 = fit0->GetErrorYhigh(i), eyl0 = fit0->GetErrorYlow(i);
        int j1 = findBin(fit1, x);    if (j1 == -1) continue;
        int jS = findBin(corrFactor, x); if (jS == -1) continue;
        double y1  = fit1->GetY()[j1], eyh1 = fit1->GetErrorYhigh(j1), eyl1 = fit1->GetErrorYlow(j1);
        double yS  = corrFactor->GetY()[jS], eyhS = corrFactor->GetErrorYhigh(jS), eylS = corrFactor->GetErrorYlow(jS);
        double ycorr = y0 + (y1 - y0) * yS;
        double ycorr_hi = ycorr, ycorr_lo = ycorr;
        for (int s0 = -1; s0 <= +1; ++s0) { for (int s1 = -1; s1 <= +1; ++s1) {
            double y0i = y0 + s0 * (s0 > 0 ? eyh0 : eyl0);
            double y1i = y1 + s1 * (s1 > 0 ? eyh1 : eyl1);
            ycorr_hi = TMath::Max(ycorr_hi, y0i + TMath::Max(0., y1i - y0i) * (yS+eyhS));
            ycorr_lo = TMath::Min(ycorr_lo, y0i + TMath::Max(0., y1i - y0i) * (yS-eyhS));
        } }

        std::cout << "x = " << x << " [" << (x-fit0->GetErrorXlow(i)) << ", " << (x+fit0->GetErrorXhigh(i)) << "] \t" <<
                     "y0  = " << y0  << " -"<<eyl0 <<"/+"<<eyh0 << " \t " <<
                     "y1  = " << y1  << " -"<<eyl1 <<"/+"<<eyh1 << " \t " <<
                     "yf  = " << yS  << " -" <<eylS<<"/+"<<eyhS<< " \t " <<
                     "yc  = " << ycorr << " -" << (ycorr-ycorr_lo) << "/+" << (ycorr_hi-ycorr) << std::endl;

        out->Set(k+1);
        out->SetPoint(k, x, ycorr);
        out->SetPointError(k, fit0->GetErrorXlow(i), fit0->GetErrorXhigh(i), ycorr - ycorr_lo, ycorr_hi - ycorr );
        k++;
    }
    return out;
}



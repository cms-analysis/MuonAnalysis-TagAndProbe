#include <TPad.h>

/** Grab a TObject from a directory
 *  even knowing only the beginning 
 *  of it's name.
 */
TObject *getFromPrefix(TDirectory *dir, TString prefixName, bool prefixOnly=true) {
    TIter next(dir->GetListOfKeys());
    TKey *key;
    if (prefixName.MaybeRegexp()) {
        TPRegexp pat(prefixName);
        while ((key = (TKey *) next())) {
            if (key->GetName() == 0) continue;
            if (TString(key->GetName()).Contains(pat)) {
                return dir->Get(key->GetName());
            }
        }
    } else {
        while ((key = (TKey *) next())) {
            if (key->GetName() == 0) continue;
            const char *match = strstr(key->GetName(), prefixName.Data());
            if (match == key->GetName() || (!prefixOnly && match != 0)) {
                return dir->Get(key->GetName());
            }
        }
    }
    return 0;
}
void plotUtil() { }

/* Other useful plotting macros */

TString preliminary = ""; //"CMS Preliminary"
TString retitle = "";
TString datalbl = "Data, xx nb^{-1}", reflbl = "Simulation";
bool autoScale = false;
bool doDiffPlot = false;
bool doRatioPlot = true;
bool doFillMC = false;
bool doPdf = true;
bool doLogX = false;
bool doSquare = false;
double yMax = 1.0;
double yMin = 0.0;
double yMaxR = 1.5;
double yMinR = 0.5;
TString extraSpam = "";

void cmsprelim() {
    TPaveText *cmsprel = new TPaveText(doSquare ? 0.40 : .55,.16,.94,.21,"NDC");
    cmsprel->SetTextSize(doSquare ? 0.040 : 0.05);
    cmsprel->SetFillColor(0);
    cmsprel->SetFillStyle(0);
    cmsprel->SetLineStyle(2);
    cmsprel->SetLineColor(0);
    cmsprel->SetTextAlign(12);
    cmsprel->SetTextFont(42);
    cmsprel->AddText(preliminary);
    cmsprel->Draw("same");
}

void doLegend(TGraphAsymmErrors *g1, TGraphAsymmErrors *g2, TString lab1, TString lab2) {
    double legend_y_offset = (preliminary != "" ? 0.07 : 0);
    double legend_y_size   = (extraSpam == "" ? .12 : .18);
    //if (g1->GetY()[g1->GetN()-1] < 0.4) {
    //    legend_y_offset = 0.75 - legend_y_size;
    //}
    TLegend *leg = new TLegend(doSquare ? .62 : .68,.15 + legend_y_offset,.92,.15 + legend_y_size + legend_y_offset);
    if (extraSpam != "") {
        leg->SetHeader(extraSpam);
        //leg->AddEntry("", extraSpam, "");
    }
    leg->AddEntry(g1, lab1, "LP");
    leg->AddEntry(g2, lab2, doFillMC ? "F" : "LP");
    leg->SetTextSize(doSquare ? 0.04 : 0.05);
    leg->SetTextFont(42);
    leg->SetShadowColor(0);
    leg->SetFillColor(0);
    if (preliminary != "") cmsprelim();
    leg->Draw();
}
void doLegend(TGraphAsymmErrors *g1, TGraphAsymmErrors *g2, TGraphAsymmErrors *g3, TString lab1, TString lab2, TString lab3) {
    double legend_y_offset = (preliminary != "" ? 0.07 : 0);
    double legend_y_size   = (extraSpam == "" ? .17 : .23);
    if (g1->GetY()[g1->GetN()-1] < 0.4) {
        legend_y_offset = 0.75 - legend_y_size;
    }
    TLegend *leg = new TLegend(doSquare ? .52 : .58,.15 + legend_y_offset,.92,.15 + legend_y_size + legend_y_offset);
    if (extraSpam != "") {
        leg->SetHeader(extraSpam);
        //leg->AddEntry("", extraSpam, "");
    }
    leg->AddEntry(g1, lab1, "LP");
    leg->AddEntry(g2, lab2, "LP");
    leg->AddEntry(g3, lab3, "LP");
    leg->SetTextSize(doSquare ? 0.04 : 0.05);
    leg->SetTextFont(42);
    leg->SetShadowColor(0);
    leg->SetFillColor(0);
    if (preliminary != "") cmsprelim();
    leg->Draw();
}


void squareCanvas(TCanvas *c) {
    c->SetCanvasSize(600,600);
    gStyle->SetPaperSize(20.,20.);
}

int findBin(TGraphAsymmErrors *g, double x) {
    for (int i = 0; i < g->GetN(); ++i) {
        double xi = g->GetX()[i];
        if ((xi - g->GetErrorXlow(i) <= x) && (x <= xi + g->GetErrorXhigh(i))) {
            return i;
        }
    }
    return -1;
}
double xmaxGraph(TGraphAsymmErrors *g) {
    double ret = 0;
    if (g != 0) {
        int n = g->GetN();
        if (n) ret =  g->GetX()[n-1] + g->GetErrorXhigh(n-1);
    }
    return ret;
}
void reTitleY(TCanvas *pl, TString ytitle) {
    TH1 *first = (TH1*) pl->GetListOfPrimitives()->At(0);
    TH1 *last = (TH1*) pl->GetListOfPrimitives()->At(pl->GetListOfPrimitives()->GetSize()-1);
    if (first) reTitleTAxis(first->GetYaxis(), ytitle);
    if (last)  reTitleTAxis(last->GetYaxis(), ytitle);
}  
void setRangeY(TCanvas *c, double min=0, double max=1.1) {
    for (size_t i = 0, n = c->GetListOfPrimitives()->GetSize(); i < n; ++i) {
        TObject *o = c->GetListOfPrimitives()->At(i);
        if (o->InheritsFrom("TH1")) ((TH1*) o)->GetYaxis()->SetRangeUser(min,max);
    }
}
const char * getXtitle(TCanvas *from) {
    TH1 *frame = (TH1*) from->GetListOfPrimitives()->At(0);
    return frame->GetXaxis()->GetTitle();
} 
void maybeLogX(TCanvas *c, TGraphAsymmErrors *h) {
    if (doLogX) {
        c->SetLogx(1);
        for (size_t i = 0, n = c->GetListOfPrimitives()->GetSize(); i < n; ++i) {
            TObject *o = c->GetListOfPrimitives()->At(i);
            if (o->InheritsFrom("TH1")) {
                TH1 *h1 = (TH1*) o;
                if (h1->GetXaxis()) {
                    h1->GetXaxis()->SetMoreLogLabels(1);
                    h1->GetXaxis()->SetTitleOffset(1.1);
                }
            }
        }
    } else {
        c->SetLogx(0);
    }
}

void reTitleTAxis(TAxis *ax, TString ytitle, double yoffset=1.3) {
   ax->SetTitle(ytitle); 
   ax->SetTitleOffset(yoffset); 
   ax->SetDecimals(true);
   ax->SetMoreLogLabels(true);
}


void doRatio(TGraphAsymmErrors *hfit, TGraphAsymmErrors *href, TString alias, const char *xtitle) {
    size_t nNZD = 0; // non-zero-denominator
    for (size_t i = 0, n = hfit->GetN(); i < n; ++i) {
        int j = findBin(href,hfit->GetX()[i]); if (j == -1) continue ;
        if (fabs(href->GetY()[j]) > 0.05) nNZD++;
    }
    TGraphAsymmErrors ratio(nNZD);
    double max = 0;
    for (size_t i = 0, k = 0, n = hfit->GetN(); i < n; ++i) {
        int j = findBin(href,hfit->GetX()[i]); if (j == -1) continue ;
        if (fabs(href->GetY()[j]) < 0.05) continue; else ++k;
        double r   = hfit->GetY()[i]/href->GetY()[j];
        double rup = (hfit->GetY()[i] == 0 ? hfit->GetErrorYhigh(i)/(href->GetY()[j]) :
                                             r*TMath::Hypot(hfit->GetErrorYhigh(i)/hfit->GetY()[i], href->GetErrorYlow(j)/href->GetY()[j]));
        double rdn = (hfit->GetY()[i] == 0 ? 0 :
                                             r*TMath::Hypot(hfit->GetErrorYlow(i)/hfit->GetY()[i],  href->GetErrorYhigh(j)/href->GetY()[j]));
        max = TMath::Max(max, fabs(r-1+rup));
        max = TMath::Max(max, fabs(r-1-rdn));
        ratio.SetPoint(k-1, hfit->GetX()[i], r);
        ratio.SetPointError(k-1, hfit->GetErrorXlow(i), hfit->GetErrorXhigh(i), rdn, rup);
    }

    ratio.Draw("AP");
    TLine line(ratio.GetX()[0]-ratio.GetErrorXlow(0), 1, ratio.GetX()[ratio.GetN()-1]+ratio.GetErrorXhigh(ratio.GetN()-1), 1);
    line.SetLineWidth(2);
    line.SetLineColor(kRed);
    line.DrawClone("SAME");
    ratio.SetLineWidth(2);
    ratio.SetLineColor(kBlack);
    ratio.SetMarkerColor(kBlack);
    ratio.SetMarkerStyle(20);
    ratio.SetMarkerSize(1.6);
    ratio.Draw("P SAME");
    if (autoScale) {
        ratio.GetYaxis()->SetRangeUser(1-1.5*max,1+1.2*max);
    } else {
        ratio.GetYaxis()->SetRangeUser(yMinR,yMaxR);
    }
    ratio.GetXaxis()->SetRangeUser(ratio.GetX()[0]-ratio.GetErrorXlow(0), ratio.GetX()[ratio.GetN()-1]+ratio.GetErrorXhigh(ratio.GetN()-1));
    ratio.GetXaxis()->SetTitle(xtitle); ratio.GetXaxis()->SetMoreLogLabels(1);
    if (datalbl) reTitleTAxis(ratio.GetYaxis(), datalbl+"/"+reflbl+" ratio");
    if (preliminary != "") cmsprelim();
    gPad->SetLogx(doLogX && (diff.GetXaxis()->GetXmin() > 0));
    gPad->Print(prefix+alias+"_ratio.png");
    if (doPdf) gPad->Print(prefix+alias+"_ratio.pdf");
}

void doDiff(RooHist *hfit, RooHist *href, TString alias, const char *xtitle) {
    double maxError = 0.7; 
    size_t nTP = 0; // non-trivial point (interval not equal to [0,1])
    for (size_t i = 0, n = hfit->GetN(); i < n; ++i) {
        int j = findBin(href,hfit->GetX()[i]); if (j == -1) continue ;
        if (hfit->GetErrorYhigh(i)+hfit->GetErrorYlow(i) >= maxError) continue;
        if (href->GetErrorYhigh(i)+href->GetErrorYlow(j) >= maxError) continue;
        nTP++;
    }
    TGraphAsymmErrors diff(nTP);
    double max = 0;
    for (size_t i = 0, n = hfit->GetN(),k=0; i < n; ++i) {
        int j = findBin(href,hfit->GetX()[i]); if (j == -1) continue ;
        if (hfit->GetErrorYhigh(i)+hfit->GetErrorYlow(i) >= maxError) continue;
        if (href->GetErrorYhigh(i)+href->GetErrorYlow(j) >= maxError) continue;
        max = TMath::Max(max, fabs(hfit->GetY()[i] - href->GetY()[j]) + fabs(hfit->GetErrorYhigh(i)) + fabs(hfit->GetErrorYlow(i)));
        max = TMath::Max(max, fabs(href->GetErrorYlow(j)) + fabs(href->GetErrorYhigh(j)));
        diff.SetPoint(k, hfit->GetX()[i], hfit->GetY()[i] - href->GetY()[j]);
        diff.SetPointError(k, hfit->GetErrorXlow(i), hfit->GetErrorXhigh(i), 
                              TMath::Hypot(hfit->GetErrorYlow(i),  href->GetErrorYlow(j)), 
                              TMath::Hypot(hfit->GetErrorYhigh(i), href->GetErrorYhigh(j))); 
        k++;
    }
    diff.SetLineWidth(2);
    diff.SetLineColor(kBlack);
    diff.SetMarkerColor(kBlack);
    diff.SetMarkerStyle(20);
    diff.SetMarkerSize(1.6);

    TLine line(diff.GetX()[0]-diff.GetErrorXlow(0), 0, diff.GetX()[diff.GetN()-1]+diff.GetErrorXhigh(diff.GetN()-1), 0);
    line.SetLineWidth(2);
    line.SetLineColor(kRed);

    diff.Draw("AP");
    line.DrawClone("SAME");
    diff.Draw("P SAME");
    diff.GetXaxis()->SetRangeUser(diff.GetX()[0]-diff.GetErrorXlow(0), diff.GetX()[diff.GetN()-1]+diff.GetErrorXhigh(diff.GetN()-1));
    if (autoScale) {
        diff.GetYaxis()->SetRangeUser(-1.5*max,1.2*max);
    } else {
        diff.GetYaxis()->SetRangeUser(-0.5,0.5);
    }
    diff.GetXaxis()->SetTitle(xtitle); //diff.GetXaxis()->SetMoreLogLabels(1);
    if (datalbl) reTitleTAxis(diff.GetYaxis(), datalbl+" - "+reflbl+" difference");
    if (preliminary != "") cmsprelim();
    //gPad->SetLogx(doLogX && (diff.GetXaxis()->GetXmin() > 0));
    gPad->Print(prefix+alias+"_diff.png");
    if (doPdf) gPad->Print(prefix+alias+"_diff.pdf");

}
/** Plot FIT from file 1 plus FIT from file 2 */
void refstack(TDirectory *fit, TDirectory *ref, TString alias, TString fitname) {
    if (ref == 0) {
        std::cerr << "REFERENCE DIR FOUND: " << fit->GetName() << std::endl;
        return;
    }

    TCanvas *pref = getFromPrefix(ref->GetDirectory("fit_eff_plots"), fitname);
    if (pref == 0) {
        std::cerr << "NOT FOUND: " << "fit_eff_plots/"+fitname << " in " << ref->GetName() << std::endl;
        return;
    }
    RooHist *href = (RooHist *) pref->FindObject("hxy_fit_eff");
    if (doFillMC) {
        href->SetLineColor(2);
        href->SetFillColor(208);
        href->SetLineStyle(0);
        href->SetMarkerColor(2);
        href->SetMarkerStyle(21);
        href->SetMarkerSize(0.4);
    } else {
        href->SetLineWidth(2);
        href->SetLineColor(kRed);
        href->SetMarkerColor(kRed);
        href->SetMarkerStyle(25);
        href->SetMarkerSize(2.0);
    }

    TCanvas *pfit = getFromPrefix(fit->GetDirectory("fit_eff_plots"), fitname);
    if (pfit == 0) {
        std::cerr << "NOT FOUND: " << "fit_eff_plots/"+fitname << " in " << fit->GetName() << std::endl;
        return;
    }
    RooHist *hfit = (RooHist *) pfit->FindObject("hxy_fit_eff");
    hfit->SetLineWidth(2);
    hfit->SetLineColor(kBlack);
    hfit->SetMarkerColor(kBlack);
    hfit->SetMarkerStyle(20);
    hfit->SetMarkerSize(1.6);

    setRangeY(pref, yMin, yMax);
    if (retitle != "") reTitleY(pref, retitle);

    pref->Draw( "" );
    if (doFillMC) href->Draw("E2 SAME");
    hfit->Draw(doFillMC ? "P0Sames" : "P SAME");
    if (datalbl) doLegend(hfit,href,datalbl,reflbl);
   
    if (doSquare) squareCanvas(pref);
    maybeLogX(pref, href); 
    gPad->Print(prefix+alias+".png");
    if (doPdf) gPad->Print(prefix+alias+".pdf");

    if (doRatioPlot) doRatio(hfit,href,alias,getXtitle(pfit)); 
    if (doDiffPlot) doDiff(hfit,href,alias,getXtitle(pfit)); 
}
/** Plot FIT from file 1 plus FIT from file 2 */
void refstackNamed(TDirectory *fit, TString alias, TString fitname, TString refname) {
    TCanvas *pref = getFromPrefix(fit->GetDirectory("fit_eff_plots"), refname);
    if (pref == 0) {
        std::cerr << "NOT FOUND: " << "fit_eff_plots/"+refname << " in " << ref->GetName() << std::endl;
        return;
    }
    RooHist *href = (RooHist *) pref->FindObject("hxy_fit_eff");
    if (doFillMC) {
        href->SetLineColor(2);
        href->SetFillColor(208);
        href->SetLineStyle(0);
        href->SetMarkerColor(2);
        href->SetMarkerStyle(21);
        href->SetMarkerSize(0.4);
    } else {
        href->SetLineWidth(2);
        href->SetLineColor(kRed);
        href->SetMarkerColor(kRed);
        href->SetMarkerStyle(25);
        href->SetMarkerSize(2.0);
    }

    TCanvas *pfit = getFromPrefix(fit->GetDirectory("fit_eff_plots"), fitname);
    if (pfit == 0) {
        std::cerr << "NOT FOUND: " << "fit_eff_plots/"+fitname << " in " << fit->GetName() << std::endl;
        return;
    }
    RooHist *hfit = (RooHist *) pfit->FindObject("hxy_fit_eff");
    hfit->SetLineWidth(2);
    hfit->SetLineColor(kBlack);
    hfit->SetMarkerColor(kBlack);
    hfit->SetMarkerStyle(20);
    hfit->SetMarkerSize(1.6);

    setRangeY(pref, yMin, yMax);
    if (retitle != "") reTitleY(pref, retitle);

    pref->Draw( "" );
    if (doFillMC) href->Draw("E2 SAME");
    hfit->Draw(doFillMC ? "P0Sames" : "P SAME");
    if (datalbl) doLegend(hfit,href,datalbl,reflbl);
   
    if (doSquare) squareCanvas(pref);
    maybeLogX(pref, href); 
    gPad->Print(prefix+alias+".png");
    if (doPdf) gPad->Print(prefix+alias+".pdf");

    if (doRatioPlot) doRatio(hfit,href,alias,getXtitle(pfit)); 
    if (doDiffPlot) doDiff(hfit,href,alias,getXtitle(pfit)); 
}

/** Plot FIT from file 1 plus FIT from file 2 plus CNT from file 3 */
void refstack3(TDirectory *fit, TDirectory *ref, TDirectory *mc, TString alias, TString fitname) {
    TCanvas *pref = getFromPrefix(ref->GetDirectory("fit_eff_plots"), fitname);
    if (pref == 0) {
        std::cerr << "NOT FOUND: " << "fit_eff_plots/"+fitname << " in " << ref->GetName() << std::endl;
        return;
    }
    RooHist *href = (RooHist *) pref->FindObject("hxy_fit_eff");
    href->SetLineWidth(2);
    href->SetLineColor(kRed);
    href->SetMarkerColor(kRed);
    href->SetMarkerStyle(25);
    href->SetMarkerSize(2.0);

    TCanvas *pmc = getFromPrefix(mc->GetDirectory("cnt_eff_plots"), fitname);
    RooHist *hmc = (RooHist *) pmc->FindObject("hxy_cnt_eff");
    hmc->SetLineWidth(2);
    hmc->SetLineColor(209);
    hmc->SetMarkerColor(209);
    hmc->SetMarkerStyle(21);
    hmc->SetMarkerSize(1.6);

    TCanvas *pfit = getFromPrefix(fit->GetDirectory("fit_eff_plots"), fitname);
    if (pfit == 0) {
        std::cerr << "NOT FOUND: " << "fit_eff_plots/"+fitname << " in " << fit->GetName() << std::endl;
        return;
    }
    RooHist *hfit = (RooHist *) pfit->FindObject("hxy_fit_eff");
    hfit->SetLineWidth(2);
    hfit->SetLineColor(kBlack);
    hfit->SetMarkerColor(kBlack);
    hfit->SetMarkerStyle(20);
    hfit->SetMarkerSize(1.6);

    setRangeY(pref, yMin, yMax);
    if (retitle != "") reTitleY(pref, retitle);

    pref->Draw("");
    hmc->Draw("P SAME");
    hfit->Draw("P SAME");
    if (datalbl) doLegend(hfit,href,hmc,"T&P "+datalbl,"T&P "+reflbl, "Simulation truth");
   
    if (doSquare) squareCanvas(pref);
    maybeLogX(pref, href); 
    gPad->Print(prefix+alias+".png");
    if (doPdf) gPad->Print(prefix+alias+".pdf");
}

/** Plot FIT from file 1 plus CNT from file 2 */
void mcstack(TDirectory *fit, TDirectory *ref, TString alias, TString name) {
    TCanvas *pref = getFromPrefix(ref->GetDirectory("cnt_eff_plots"), name);
    RooHist *href = (RooHist *) pref->FindObject("hxy_cnt_eff");
    href->SetLineWidth(2);
    href->SetLineColor(209);
    href->SetMarkerColor(209);
    href->SetMarkerStyle(25);
    href->SetMarkerSize(2.0);

    TCanvas *pfit = getFromPrefix(fit->GetDirectory("fit_eff_plots"), name);
    RooHist *hfit = (RooHist *) pfit->FindObject("hxy_fit_eff");
    hfit->SetLineWidth(2);
    hfit->SetLineColor(206);
    hfit->SetMarkerColor(206);
    hfit->SetMarkerStyle(20);
    hfit->SetMarkerSize(1.6);

    setRangeY(pref, yMin, yMax);
    if (doSquare) squareCanvas(pref);
    if (retitle != "") reTitleY(pref, retitle);
    pref->Draw();
    hfit->Draw("P SAME");

    doLegend(hfit,href,datalbl,reflbl);
    maybeLogX(pref, href); 
    gPad->Print(prefix+alias+".png");
    if (doPdf) gPad->Print(prefix+alias+".pdf");

    if (doRatioPlot) doRatio(hfit,href,alias,getXtitle(pfit)); 
    if (doDiffPlot) doDiff(hfit,href,alias,getXtitle(pfit)); 
}


TGraphAsymmErrors *getFit(TDirectory *fit, TString fitname) {
    TCanvas *pfit = getFromPrefix(fit->GetDirectory("fit_eff_plots"), fitname);
    if (pfit == 0) {
        std::cerr << "NOT FOUND: " << "fit_eff_plots/"+fitname << " in " << fit->GetName() << std::endl;
        return;
    }
    if (retitle != "") reTitleY(pfit, retitle);

    RooHist *hfit = (RooHist *) pfit->FindObject("hxy_fit_eff");
    if (hfit == 0) {
        std::cerr << "NOT FOUND: " << "fit_eff_plots/"+fitname << "/hxy_fit_eff in " << fit->GetName() << std::endl;
        pfit->ls();
        return;
    }
    return hfit;
}
/** Plot just one set */
TGraphAsymmErrors *single( TDirectory *fit, TString alias, TString fitname) {
    TCanvas *pfit = getFromPrefix(fit->GetDirectory("fit_eff_plots"), fitname);
    if (pfit == 0) {
        std::cerr << "NOT FOUND: " << "fit_eff_plots/"+fitname << " in " << fit->GetName() << std::endl;
        return;
    }
    if (retitle != "") reTitleY(pfit, retitle);

    RooHist *hfit = (RooHist *) pfit->FindObject("hxy_fit_eff");
    if (hfit == 0) {
        std::cerr << "NOT FOUND: " << "fit_eff_plots/"+fitname << "/hxy_fit_eff in " << fit->GetName() << std::endl;
        pfit->ls();
        return;
    }
    hfit->SetLineWidth(2);
    hfit->SetLineColor(kBlue);
    hfit->SetMarkerColor(kBlue);
    hfit->SetMarkerStyle(20);
    hfit->SetMarkerSize(1.6);

    setRangeY(pfit, yMin, yMax);
    if (doSquare) squareCanvas(pfit);
    if (preliminary != "") cmsprelim();
    maybeLogX(pfit, hfit); 
    pfit->Print(prefix+alias+".png"); 
    if (doPdf) pfit->Print(prefix+alias+".pdf"); 

    return (TGraphAsymmErrors *) hfit->Clone();
}

void EffPalette()
{
   static Int_t  colors[50];
   static Bool_t initialized = kFALSE;

   double mid = 0.5*(yMin+yMax);
   Double_t Red[3]    = { 1.00, 1.00, 0.00 };
   Double_t Green[3]  = { 0.00, 1.00, 1.00 };
   Double_t Blue[3]   = { 0.00, 0.00, 0.00 };
   Double_t Length[3] = { 0,    0.75,  1.0  };

   if(!initialized){
      Int_t FI = TColor::CreateGradientColorTable(3,Length,Red,Green,Blue,50);
      for (int i=0; i<50; i++) colors[i] = FI+i;
      initialized = kTRUE;
   }
   //gStyle->SetPalette(50,colors); // not needed??
}


/** Plot just one set */
void single2D( TDirectory *fit, TString alias, TString fitname) {
    gStyle->SetPaintTextFormat(".4f");
    TCanvas *pfit = getFromPrefix(fit->GetDirectory("fit_eff_plots"), fitname);
    if (pfit == 0) {
        std::cerr << "NOT FOUND: " << "fit_eff_plots/"+fitname << " in " << fit->GetName() << std::endl;
        return;
    }
    TH2 *hfit = pfit->FindObject(pfit->GetName());
    hfit->SetTitle(retitle == "" ?  "efficiency" : retitle);
    hfit->GetZaxis()->SetTitle("");
    hfit->GetZaxis()->SetRangeUser(yMin, yMax);

    EffPalette();
    //gStyle->SetPalette(99);

    c1->cd();
    double orm = c1->GetRightMargin();
    c1->SetRightMargin(0.14);
    hfit->Draw();
    if (doSquare) squareCanvas(c1);
    if (preliminary != "") cmsprelim();
    pfit->Print(prefix+alias+".png"); 
    if (doPdf) pfit->Print(prefix+alias+".pdf"); 
    c1->SetRightMargin(orm);
}

/** Reset line styles and colors, which get messed up by tdrStyle */
void prettyLine(TCanvas *canv, int pad, const char *cname, int color) {
    RooCurve *c = (RooCurve *) canv->GetPad(pad)->FindObject(cname);
    c->SetLineWidth(2);
    c->SetLineColor(color);
}
void prettyLines(TCanvas *c) {
   prettyLine(c, 1, "pdfPass_Norm[mass]",                      kRed  );
   prettyLine(c, 1, "pdfPass_Norm[mass]_Comp[backgroundPass]", kBlue );
   prettyLine(c, 2, "pdfFail_Norm[mass]",                      kRed  );
   prettyLine(c, 2, "pdfFail_Norm[mass]_Comp[backgroundFail]", kBlue );
   prettyLine(c, 3, "simPdf_Norm[mass]",                                     kRed  );
   prettyLine(c, 3, "simPdf_Norm[mass]_Comp[backgroundPass,backgroundFail]", kBlue );
}
void doCanvas(TDirectory *dir, int binsx, int binsy, const char * easyname, const char * truename) {
    gSystem->mkdir(prefix+"canvases/",true);
    char buff[1023], baff[1023];
    for (int i = 0; i < binsx; ++i) {
        for (int j = 0; j < binsy; ++j) {
            if (binsx != 1 && binsy != 1) {
                sprintf(buff,easyname,i,j);
                sprintf(baff,truename,i,j);
            } else if (binsx != 1) {
                sprintf(buff,easyname,i);
                sprintf(baff,truename,i);
            } else if (binsy != 1) {
                sprintf(buff,easyname,j);
                sprintf(baff,truename,j);
            } else {
                sprintf(buff,easyname);
                sprintf(baff,truename);
            }
            TDirectory *subdir = (TDirectory *) getFromPrefix(dir, baff, /*prefixOnly=*/false);
            if (subdir == 0) {
                //std::cerr << "Didn't find '" << baff << "*' in " << dir->GetName() << std::endl;
                continue;
            }
            TCanvas *fitc = (TCanvas *) subdir->Get("fit_canvas");
            if (fitc == 0) {
                std::cerr << "Didn't find " << TString(baff) << "/fit_canvas in " << dir->GetName() << std::endl;
                continue;
            }
            fitc->Draw(); 
            prettyLines(fitc);
            fitc->Print(prefix+TString("canvases/")+buff+"_fit.png");
        }
    }
}

enum GraphMergeAlgo { Combine, LowestUnc, Last };
TGraphAsymmErrors *mergeGraphs(TGraphAsymmErrors *g1, TGraphAsymmErrors *g2, TGraphAsymmErrors *g3=0, TGraphAsymmErrors *g4=0, GraphMergeAlgo algo=LowestUnc, double yMin = 0.001) {
    TGraphAsymmErrors *graphs[4];
    graphs[0] = g1;
    graphs[1] = g2;
    graphs[2] = g3;
    graphs[3] = g4;
    int ng = 0; 
    for (int i = 0; i < 4; ++i) { if (graphs[i]) ng++; else break; }
    int n = (g1 ? g1->GetN() : 0);
    TGraphAsymmErrors *ret = new TGraphAsymmErrors(n);
    for (int i = 0; i < n; ++i) {
        double syw2 = 0, sw2 = 0, slo = 0, shi = 0;
        double xbin = g1->GetX()[i], sxw2 = 0;
        double xmax = xbin + g1->GetErrorXhigh(i), xmin = xbin - g1->GetErrorXlow(i);
        double ytry = 0, ntry = 0;
        for (int j = 0; j < ng; ++j) {
            int k = (j == 0 ? i : findBin(graphs[j], xbin));
            if (k == -1) continue;
            double y = graphs[j]->GetY()[k];
            ytry += y; ntry += 1;
        }
        for (int j = 0; j < ng; ++j) {
            int k = (j == 0 ? i : findBin(graphs[j], xbin));
            if (k == -1) continue;
            double x = graphs[j]->GetX()[k];
            double y = graphs[j]->GetY()[k];
            if (algo == Combine && ytry/ntry > yMin && y < yMin) continue;
            double ylo = graphs[j]->GetErrorYlow(k);
            double yhi = graphs[j]->GetErrorYhigh(k);
            double w2 = 1.0/TMath::Max(ylo,yhi); w2 *= w2;
            //std::cout << "For graph " << j << " use point  " << k << ", x = " << x << ", y = " << y << "  -" << ylo << "/+" << yhi  << std::endl;
            switch (algo) {
                case Combine:
                    sw2  += w2;
                    syw2 += w2*y;
                    sxw2 += w2*x;
                    shi  += yhi*yhi*w2*w2;
                    slo  += ylo*ylo*w2*w2;
                    break;
                case LowestUnc:
                    if (w2 > sw2) {
                        sw2  = w2;
                        syw2 = w2*y;
                        sxw2 = w2*x;
                        shi  = yhi*yhi*w2*w2;
                        slo  = ylo*ylo*w2*w2;
                    }
                    break;
                case Last:
                    sw2  = w2;
                    syw2 = w2*y;
                    sxw2 = w2*x;
                    shi  = yhi*yhi*w2*w2;
                    slo  = ylo*ylo*w2*w2;
                    break;
            }
        }
        ret->SetPoint(i, sxw2/sw2, syw2/sw2);
        ret->SetPointError(i, sxw2/sw2 - xmin, xmax - sxw2/sw2, sqrt(shi)/sw2, sqrt(slo)/sw2);
    }
    return ret;
}    



#include <TCanvas.h>
#include <TColor.h>
#include <TFile.h>
#include <TGraphAsymmErrors.h>
#include <TH1.h>
#include <TH2.h>
#include <TKey.h>
#include <TLegend.h>
#include <TMath.h>
#include <TPRegexp.h>
#include <TPad.h>
#include <TPaveText.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TStyle.h>
#include <TSystem.h>
#include <math.h>
#include <stdlib.h>
#include <iostream>

TString prefix = "plots_dev/muonid/";
TString basedir  = "tpTree";

TFile *ref = 0;

TCanvas *c1 = 0;


/** Grab a TObject from a directory
 *  even knowing only the beginning 
 *  of it's name.
 */
TObject *getFromPrefix(TDirectory *dir, TString prefixName, bool prefixOnly=true) {
    if (dir == 0) { std::cerr << "ERROR: trying to get " << prefixName << " from null dir." << std::endl; return 0; }
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
TString retitle = "", retitleX = "";
TString datalbl = "Data, xx nb^{-1}", reflbl = "Simulation";
bool autoScale = false;
bool doDiffPlot = false;
bool doRatioPlot = true;
bool doFillMC = false;
bool doPdf = true;
bool doEps = true;
bool doTxt = true;
bool doLogX = false;
bool doSquare = false;
bool doSanity = false;
double yMax = 1.0;
double yMin = 0.0;
double yMaxR = 1.5;
double yMinR = 0.5;
double yMaxD = +0.5;
double yMinD = -0.5;
TString extraSpam = "";
TFile *fOut = 0;

int lineColorRef_FillMC = 62;
int fillColorRef_FillMC = 65;

double cmsprel_xoffs = 0;
double cmsprel_yoffs = 0;
void cmsprelim(double xoffs=-99,double yoffs=-99) {
    if (xoffs == -99) xoffs = cmsprel_xoffs;
    if (yoffs == -99) yoffs = cmsprel_yoffs;
    TPaveText *cmsprel = new TPaveText(xoffs+(doSquare ? 0.40 : .55),yoffs+.16,xoffs+.94,yoffs+.21,"NDC");
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

double doLegend_xoffs = 0;
double doLegend_yoffs = 0;
void doLegend(TGraphAsymmErrors *g1, TGraphAsymmErrors *g2, TString lab1, TString lab2) {
    double legend_y_offset = doLegend_yoffs+(preliminary != "" ? 0.07 : 0);
    double legend_y_size   = (extraSpam == "" ? .12 : .18);
    //if (g1->GetY()[g1->GetN()-1] < 0.4) {
    //    legend_y_offset = 0.75 - legend_y_size;
    //}
    TLegend *leg = new TLegend(doLegend_xoffs+(doSquare ? .62 : .68),.15 + legend_y_offset,doLegend_xoffs+.92,.15 + legend_y_size + legend_y_offset);
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
    double legend_y_offset = doLegend_yoffs+(preliminary != "" ? 0.07 : 0);
    double legend_y_size   = (extraSpam == "" ? .17 : .23);
    //if (g1->GetY()[g1->GetN()-1] < 0.4) {
    //    legend_y_offset = 0.75 - legend_y_size;
    //}
    TLegend *leg = new TLegend(doLegend_xoffs+(doSquare ? .52 : .58),.15 + legend_y_offset,doLegend_xoffs+.92,.15 + legend_y_size + legend_y_offset);
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
void doLegend(TGraphAsymmErrors *g1, TGraphAsymmErrors *g2, TGraphAsymmErrors *g3, TGraphAsymmErrors *g4, TString lab1, TString lab2, TString lab3, TString lab4) {
    double legend_y_offset = (preliminary != "" ? 0.07 : 0);
    double legend_y_size   = (extraSpam == "" ? .22 : .28);
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
    leg->AddEntry(g4, lab4, "LP");
    leg->SetTextSize(doSquare ? 0.04 : 0.05);
    leg->SetTextFont(42);
    leg->SetShadowColor(0);
    leg->SetFillColor(0);
    if (preliminary != "") cmsprelim();
    leg->Draw();
}


void squareCanvas(TVirtualPad *c) {
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
double xminGraph(TGraphAsymmErrors *g) {
    double ret = 0;
    if (g != 0) {
        int n = g->GetN();
        if (n) ret =  g->GetX()[0] - g->GetErrorXlow(0);
    }
    return ret;
}

double xmaxGraph(TGraphAsymmErrors *g) {
    double ret = 0;
    if (g != 0) {
        int n = g->GetN();
        if (n) ret =  g->GetX()[n-1] + g->GetErrorXhigh(n-1);
    }
    return ret;
}
double ymaxGraph(TGraphAsymmErrors *g) {
    double ret = 0;
    if (g != 0) {
        int n = g->GetN();
        for (int i = 0; i < n; ++i) {
            ret = TMath::Max(ret, g->GetY()[i] + g->GetErrorYhigh(i));
        }
    }
    return ret;
}

void reTitleTAxis(TAxis *ax, TString ytitle, double yoffset=1.3) {
   ax->SetTitle(ytitle); 
   ax->SetTitleOffset(yoffset); 
   ax->SetDecimals(true);
   ax->SetMoreLogLabels(true);
}

void reTitleY(TCanvas *pl, TString ytitle) {
    TH1 *first = (TH1*) pl->GetListOfPrimitives()->At(0);
    //TH1 *last = (TH1*) pl->GetListOfPrimitives()->At(pl->GetListOfPrimitives()->GetSize()-1);
    if (first) reTitleTAxis(first->GetYaxis(), ytitle);
    //if (last)  reTitleTAxis(last->GetYaxis(), ytitle);
}  
void reTitleX(TCanvas *pl, TString xtitle) {
    TH1 *first = (TH1*) pl->GetListOfPrimitives()->At(0);
    //TH1 *last = (TH1*) pl->GetListOfPrimitives()->At(pl->GetListOfPrimitives()->GetSize()-1);
    if (first) reTitleTAxis(first->GetXaxis(), xtitle, 0.9);
    //if (last)  reTitleTAxis(last->GetXaxis(), xtitle, 0.9);
}  
void drawFrame(TCanvas *pl) {
    TH1 *first = (TH1*) pl->GetListOfPrimitives()->At(0);
    //TH1 *last = (TH1*) pl->GetListOfPrimitives()->At(pl->GetListOfPrimitives()->GetSize()-1);
    first->Draw("");
    //last->Draw("SAME");
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

void sanityCheck(TGraphAsymmErrors *fit) {
    for (int i = 0; i < fit->GetN(); ++i) {
        double x = fit->GetX()[i], y = fit->GetY()[i], eyl = fit->GetErrorYlow(i), eyh = fit->GetErrorYhigh(i);
        if (eyl == 0 && y > 0) {
            fprintf(stdout,"Error: %s, x = %g: y = %g  -%g/+%g\n", fit->GetName(), x,y,eyl,eyh);
            if (y < 0.5 && (y-3*eyl < 0)) { eyl = y; } else { eyl = eyh; }
        }
        if (eyl == 0 && y == 1) {
            fprintf(stdout,"Error: %s, x = %g: y = %g  -%g/+%g\n", fit->GetName(), x,y,eyl,eyh);
            if (y > 0.5 && (y + 3*eyh > 0)) { eyh = 1-y; } else { eyh = eyl; }
        }
        if ((eyl > 3*eyh) && (y+eyh < 0.9997 || eyh > 0.01)) {
            fprintf(stdout,"Error: %s, x = %g: y = %g  -%g/+%g: big eyl\n", fit->GetName(), x,y,eyl,eyh);
            eyl = eyh;
            
        }
        if ((eyh > 3*eyl) && (y-eyl > 0.0003 || eyl > 0.01)) {
            fprintf(stdout,"Error: %s, x = %g: y = %g  -%g/+%g: big eyh\n", fit->GetName(), x,y,eyl,eyh);
            eyh = eyl;
        }
        if (eyh > 0.2) eyh = TMath::Min(1-y, (i > 0 ? fit->GetErrorYhigh(i-1) : 0.01));
        if (eyl > 0.2) eyl = TMath::Min(  y, (i > 0 ? fit->GetErrorYlow(i-1)  : 0.01));
        if (i > 0 && i < fit->GetN()-1) {
            double y1 = fit->GetY()[i-1];
            double y2 = fit->GetY()[i+1];
            if (fabs(y - y1) > 0.2 && fabs(y - y2) > 0.2 && ((y - y1)*(y - y2) > 0)) {
                fit->SetPoint(i, x, 0.5*(y1+y2));
            }
        }
        fit->SetPointEYhigh(i, eyh);
        fit->SetPointEYlow(i, eyl);
    }
    fflush(stdout);
}

void printGraph(TGraphAsymmErrors *graph, TString alias) {
    if (graph == 0 || graph->GetN() == 0 || !graph->InheritsFrom("TGraphAsymmErrors")) return;
    FILE *text = fopen((prefix+alias+".txt").Data(), "w");
    if (text == 0) return;
    fprintf(text, " %7s  %7s  %7s    %6s   -%4s/ +%4s\n","  x min","  x avg","  x max","  eff ","err ","err ");
    fprintf(text, " %7s  %7s  %7s    %6s  %6s/%6s\n"," ------"," ------"," ------"," -----","------","------");
    for (int i = 0; i < graph->GetN(); ++i) {
        double x = graph->GetX()[i], xlo = graph->GetErrorXlow(i), xhi = graph->GetErrorXhigh(i);
        double y = graph->GetY()[i], ylo = graph->GetErrorYlow(i), yhi = graph->GetErrorYhigh(i);
        fprintf(text, " %7.3f  %7.3f  %7.3f    %6.2f  -%5.2f/+%5.2f\n", x-xlo,x,x+xhi,y*100,ylo*100,yhi*100);
    }
    fclose(text);
}
void printGraphs(TGraphAsymmErrors *graph, TGraphAsymmErrors *graph2, TString alias) {
    if (graph == 0 || graph->GetN() == 0) return;
    FILE *text = fopen((prefix+alias+".txt").Data(), "w");
    fprintf(text, " %7s  %7s  %7s    %6s   -%4s/ +%4s    %6s   -%4s/ +%4s\n","  x min","  x avg","  x max","  eff "," err"," err","  ref "," err"," err");
    fprintf(text, " %7s  %7s  %7s    %6s  %6s/%6s    %6s  %6s/%6s\n"," ------"," ------"," ------"," -----","------","------"," -----","------","------");
    for (int i = 0; i < graph->GetN(); ++i) {
        double x = graph->GetX()[i], xlo = graph->GetErrorXlow(i), xhi = graph->GetErrorXhigh(i);
        double y = graph->GetY()[i], ylo = graph->GetErrorYlow(i), yhi = graph->GetErrorYhigh(i);
        int j = findBin(graph2, x);
        if (j == -1) {
            fprintf(text, " %7.3f  %7.3f  %7.3f    %6.2f  -%5.2f/+%5.2f    %6s  %6s/%6s\n", x-xlo,x,x+xhi,y*100,ylo*100,yhi*100,"--.--","-.--","-.--");
        } else {
            double y2 = graph2->GetY()[i], y2lo = graph2->GetErrorYlow(i), y2hi = graph2->GetErrorYhigh(i);
            fprintf(text, " %7.3f  %7.3f  %7.3f    %6.2f  -%5.2f/+%5.2f    %6.2f  -%5.2f/+%5.2f\n", x-xlo,x,x+xhi,y*100,ylo*100,yhi*100, y2*100,y2lo*100,y2hi*100);
        }
    }
    fclose(text);
}



void doRatio(TGraphAsymmErrors *hfit, TGraphAsymmErrors *href, TString alias, const char *xtitle) {
    if (hfit->GetN() == 0 || href->GetN() == 0) return;
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

    TH1F *frame = new TH1F("frameratio","frameratio",1, ratio.GetX()[0]-ratio.GetErrorXlow(0), ratio.GetX()[ratio.GetN()-1]+ratio.GetErrorXhigh(ratio.GetN()-1));
    frame->Draw("");
    frame->GetXaxis()->SetTitle(xtitle); frame->GetXaxis()->SetMoreLogLabels(1); frame->GetXaxis()->SetNoExponent(1);
    TLine line(ratio.GetX()[0]-ratio.GetErrorXlow(0), 1, ratio.GetX()[ratio.GetN()-1]+ratio.GetErrorXhigh(ratio.GetN()-1), 1);
    line.SetLineWidth(2);
    line.SetLineColor(fillColorRef_FillMC);
    line.DrawClone("SAME");
    ratio.SetLineWidth(3);
    ratio.SetLineColor(kBlack);
    ratio.SetMarkerColor(kBlack);
    ratio.SetMarkerStyle(20);
    ratio.SetMarkerSize(1.6);
    ratio.Draw("P SAME");
    if (autoScale) {
        frame->GetYaxis()->SetRangeUser(1-1.5*max,1+1.2*max);
    } else {
        frame->GetYaxis()->SetRangeUser(yMinR,yMaxR);
    }
    if (datalbl) reTitleTAxis(frame->GetYaxis(), datalbl+" / "+reflbl+" ratio");
    if (preliminary != "") cmsprelim();
    if (doSquare) squareCanvas(gPad);
    if (yMaxR - yMinR < 0.05) frame->GetYaxis()->SetTitleOffset(1.6);
    gPad->SetLogx(doLogX && (frame->GetXaxis()->GetXmin() > 0));
    gPad->Print(prefix+alias+"_ratio.png");
    if (doPdf) gPad->Print(prefix+alias+"_ratio.pdf");
    if (doTxt) printGraph(&ratio,alias+"_ratio");
    delete frame;
}

void doDiff(TGraphAsymmErrors *hfit, TGraphAsymmErrors *href, TString alias, const char *xtitle) {
    if (hfit->GetN() == 0 || href->GetN() == 0) return;
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

    TH1F *frame = new TH1F("frameDiff","frameDiff",1, diff.GetX()[0]-diff.GetErrorXlow(0), diff.GetX()[diff.GetN()-1]+diff.GetErrorXhigh(diff.GetN()-1));
    frame->Draw("");
    frame->GetXaxis()->SetTitle(xtitle); frame->GetXaxis()->SetMoreLogLabels(1); frame->GetXaxis()->SetNoExponent(1);
    if (autoScale) {
        frame->GetYaxis()->SetRangeUser(-1.5*max,1.2*max);
    } else {
        frame->GetYaxis()->SetRangeUser(yMinD, yMaxD);
    }
    line.DrawClone("SAME");
    diff.Draw("P SAME");
    frame->GetXaxis()->SetTitle(xtitle); //diff.GetXaxis()->SetMoreLogLabels(1);
    if (datalbl) reTitleTAxis(frame->GetYaxis(), datalbl+" - "+reflbl+" difference");
    if (preliminary != "") cmsprelim();
    if (doSquare) squareCanvas(gPad);
    gPad->SetLogx(doLogX && (frame->GetXaxis()->GetXmin() > 0));
    gPad->Print(prefix+alias+"_diff.png");
    if (doPdf) gPad->Print(prefix+alias+"_diff.pdf");
    if (doTxt) printGraph(&diff,alias+"_diff");
    delete frame;
}
/** Plot FIT from file 1 plus FIT from file 2 */
void refstack(TDirectory *fit, TDirectory *refd, TString alias, TString fitname) {
    if (fit == 0) {
        std::cerr << "ERROR: refstack called with missing dirs: alias = " << alias << std::endl;
    }
    if (refd == 0) {
        std::cerr << "REFERENCE NOT DIR FOUND FOR: " << fit->GetName() << std::endl;
        return;
    }

    TCanvas *pref = (TCanvas*) getFromPrefix(refd->GetDirectory("fit_eff_plots"), fitname);
    if (pref == 0) {
        std::cerr << "NOT FOUND: " << "fit_eff_plots/"+fitname << " in " << refd->GetName() << std::endl;
        return;
    }
    TGraphAsymmErrors *href = (TGraphAsymmErrors *) pref->FindObject("hxy_fit_eff");
    if (doFillMC) {
        href->SetLineColor(lineColorRef_FillMC);
        href->SetFillColor(fillColorRef_FillMC);
        href->SetLineStyle(0);
        href->SetMarkerColor(lineColorRef_FillMC);
        href->SetMarkerStyle(21);
        href->SetMarkerSize(0.4);
    } else {
        href->SetLineWidth(2);
        href->SetLineColor(kRed);
        href->SetMarkerColor(kRed);
        href->SetMarkerStyle(25);
        href->SetMarkerSize(2.0);
    }

    TCanvas *pfit = (TCanvas*) getFromPrefix(fit->GetDirectory("fit_eff_plots"), fitname);
    if (pfit == 0) {
        std::cerr << "NOT FOUND: " << "fit_eff_plots/"+fitname << " in " << fit->GetName() << std::endl;
        return;
    }
    TGraphAsymmErrors *hfit = (TGraphAsymmErrors *) pfit->FindObject("hxy_fit_eff");
    hfit->SetLineWidth(2);
    hfit->SetLineColor(kBlack);
    hfit->SetMarkerColor(kBlack);
    hfit->SetMarkerStyle(20);
    hfit->SetMarkerSize(1.6);

    setRangeY(pref, yMin, yMax);
    if (retitle  != "") reTitleY(pref, retitle);
    if (retitleX != "") reTitleX(pref, retitleX);

    if (doSanity) sanityCheck(href);
    if (doSanity) sanityCheck(hfit);

    if (doFillMC) {
        if (c1 == 0) c1 = new TCanvas("c1","c1");
        c1->cd(); c1->Clear();
        drawFrame(pref);
        href->Draw("E2 SAME");
    } else {
        pref->Draw( "" );
    }
    hfit->Draw(doFillMC ? "P0Sames" : "P SAME");
    if (datalbl) doLegend(hfit,href,datalbl,reflbl);
   
    if (doSquare) squareCanvas(pref);
    maybeLogX(pref, href); 
    gPad->Print(prefix+alias+".png");
    if (doPdf) gPad->Print(prefix+alias+".pdf");
    if (doEps) gPad->Print(prefix+alias+".eps");

    if (fOut) { fOut->WriteTObject(hfit,"fit_"+alias, "Overwrite"); fOut->WriteTObject(href,"ref_"+alias, "Overwrite"); }
    if (doTxt)  printGraph(hfit,"fit_"+alias);
    if (doTxt)  printGraph(href,"ref_"+alias);
    if (doTxt)  printGraphs(hfit, href, alias);
    // move these below the doTxt
    if (doRatioPlot) doRatio(hfit,href,alias,getXtitle(pfit)); 
    if (doDiffPlot) doDiff(hfit,href,alias,getXtitle(pfit)); 
}
/** Plot FIT from file 1 plus FIT from file 2 */
void refstackNamed(TDirectory *fit, TString alias, TString fitname, TString refname) {
    TCanvas *pref = (TCanvas*) getFromPrefix(fit->GetDirectory("fit_eff_plots"), refname);
    if (pref == 0) {
        std::cerr << "NOT FOUND: " << "fit_eff_plots/"+refname << " in " << fit->GetName() << std::endl;
        return;
    }
    TGraphAsymmErrors *href = (TGraphAsymmErrors *) pref->FindObject("hxy_fit_eff");
    if (doFillMC) {
        href->SetLineColor(lineColorRef_FillMC);
        href->SetFillColor(fillColorRef_FillMC);
        href->SetLineStyle(0);
        href->SetMarkerColor(lineColorRef_FillMC);
        href->SetMarkerStyle(21);
        href->SetMarkerSize(0.4);
    } else {
        href->SetLineWidth(2);
        href->SetLineColor(kRed);
        href->SetMarkerColor(kRed);
        href->SetMarkerStyle(25);
        href->SetMarkerSize(2.0);
    }

    TCanvas *pfit = (TCanvas*) getFromPrefix(fit->GetDirectory("fit_eff_plots"), fitname);
    if (pfit == 0) {
        std::cerr << "NOT FOUND: " << "fit_eff_plots/"+fitname << " in " << fit->GetName() << std::endl;
        return;
    }
    TGraphAsymmErrors *hfit = (TGraphAsymmErrors *) pfit->FindObject("hxy_fit_eff");
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
    if (doEps) gPad->Print(prefix+alias+".eps");

    if (doRatioPlot) doRatio(hfit,href,alias,getXtitle(pfit)); 
    if (doDiffPlot) doDiff(hfit,href,alias,getXtitle(pfit)); 
    if (fOut) { fOut->WriteTObject(hfit,"fit_"+alias, "Overwrite"); fOut->WriteTObject(href,"ref_"+alias, "Overwrite"); }
}

/** Plot FIT from file 1 plus FIT from file 2 plus CNT from file 3 */
void refstack3(TDirectory *fit, TDirectory *refd, TDirectory *mc, TString alias, TString fitname) {
    TCanvas *pref = (TCanvas *) getFromPrefix(refd->GetDirectory("fit_eff_plots"), fitname);
    if (pref == 0) {
        std::cerr << "NOT FOUND: " << "fit_eff_plots/"+fitname << " in " << refd->GetName() << std::endl;
        return;
    }
    TGraphAsymmErrors *href = (TGraphAsymmErrors *) pref->FindObject("hxy_fit_eff");
    href->SetLineWidth(2);
    href->SetLineColor(kRed);
    href->SetMarkerColor(kRed);
    href->SetMarkerStyle(25);
    href->SetMarkerSize(2.0);

    TCanvas *pmc = (TCanvas *) getFromPrefix(mc->GetDirectory("cnt_eff_plots"), fitname);
    TGraphAsymmErrors *hmc = (TGraphAsymmErrors *) pmc->FindObject("hxy_cnt_eff");
    hmc->SetLineWidth(2);
    hmc->SetLineColor(209);
    hmc->SetMarkerColor(209);
    hmc->SetMarkerStyle(21);
    hmc->SetMarkerSize(1.6);

    TCanvas *pfit = (TCanvas *) getFromPrefix(fit->GetDirectory("fit_eff_plots"), fitname);
    if (pfit == 0) {
        std::cerr << "NOT FOUND: " << "fit_eff_plots/"+fitname << " in " << fit->GetName() << std::endl;
        return;
    }
    TGraphAsymmErrors *hfit = (TGraphAsymmErrors *) pfit->FindObject("hxy_fit_eff");
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
    if (doEps) gPad->Print(prefix+alias+".eps");
    if (fOut) { fOut->WriteTObject(hfit,"fit_"+alias, "Overwrite"); fOut->WriteTObject(href,"ref_"+alias, "Overwrite"); fOut->WriteTObject(hmc,"mct_"+alias, "Overwrite"); }
}

/** Plot FIT from file 1 plus CNT from file 2 */
void mcstack(TDirectory *fit, TDirectory *refd, TString alias, TString name) {
    TCanvas *pref = (TCanvas *) getFromPrefix(refd->GetDirectory("cnt_eff_plots"), name);
    TGraphAsymmErrors *href = (TGraphAsymmErrors *) pref->FindObject("hxy_cnt_eff");
    href->SetLineWidth(2);
    href->SetLineColor(209);
    href->SetMarkerColor(209);
    href->SetMarkerStyle(25);
    href->SetMarkerSize(2.0);

    TCanvas *pfit = (TCanvas *) getFromPrefix(fit->GetDirectory("fit_eff_plots"), name);
    TGraphAsymmErrors *hfit = (TGraphAsymmErrors *) pfit->FindObject("hxy_fit_eff");
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
    if (doEps) gPad->Print(prefix+alias+".eps");

    if (doRatioPlot) doRatio(hfit,href,alias,getXtitle(pfit)); 
    if (doDiffPlot) doDiff(hfit,href,alias,getXtitle(pfit)); 
    if (fOut) { fOut->WriteTObject(hfit,"fit_"+alias, "Overwrite"); fOut->WriteTObject(href,"mc_"+alias, "Overwrite"); }
}


TGraphAsymmErrors *getFit(TDirectory *fit, TString fitname, TString algo="fit") {
    if (fit == 0) { std::cerr << "getFit called with " << fitname << " and fit = 0" << std::endl; return 0; }
    TCanvas *pfit = (TCanvas *) getFromPrefix(fit->GetDirectory(algo+"_eff_plots"), fitname);
    if (pfit == 0) {
        std::cerr << "NOT FOUND: " << algo+"_eff_plots/"+fitname << " in " << fit->GetName() << std::endl;
        return 0;
    }
    if (retitle != "") reTitleY(pfit, retitle);

    TGraphAsymmErrors *hfit = (TGraphAsymmErrors *) pfit->FindObject("hxy_"+algo+"_eff");
    if (hfit == 0) {
        std::cerr << "NOT FOUND: " << algo+"_eff_plots/"+fitname << "/hxy_"+algo+"_eff in " << fit->GetName() << std::endl;
        pfit->ls();
        return 0;
    }
    return hfit;
}
/** Plot just one set */
TGraphAsymmErrors *single( TDirectory *fit, TString alias, TString fitname) {
    if (fit == 0) {
        std::cerr << "ERROR: signle called with missing dirs: alias = " << alias << std::endl;
    }
    TCanvas *pfit = (TCanvas *) getFromPrefix(fit->GetDirectory("fit_eff_plots"), fitname);
    if (pfit == 0) {
        std::cerr << "NOT FOUND: " << "fit_eff_plots/"+fitname << " in " << fit->GetName() << std::endl;
        return 0;
    }
    if (retitle != "") reTitleY(pfit, retitle);

    TGraphAsymmErrors *hfit = (TGraphAsymmErrors *) pfit->FindObject("hxy_fit_eff");
    if (hfit == 0) {
        std::cerr << "NOT FOUND: " << "fit_eff_plots/"+fitname << "/hxy_fit_eff in " << fit->GetName() << std::endl;
        pfit->ls();
        return 0;
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
    if (doEps) pfit->Print(prefix+alias+".eps");

    if (fOut) { fOut->WriteTObject(hfit,"fit_"+alias, "Overwrite"); }
    if (doTxt)  printGraph(hfit,alias);
    return (TGraphAsymmErrors *) hfit->Clone();
}

void EffPalette()
{
   static Int_t  colors[50];
   static Bool_t initialized = kFALSE;

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
    TCanvas *pfit = (TCanvas*) getFromPrefix(fit->GetDirectory("fit_eff_plots"), fitname);
    if (pfit == 0) {
        std::cerr << "NOT FOUND: " << "fit_eff_plots/"+fitname << " in " << fit->GetName() << std::endl;
        return;
    }
    TH2 *hfit = (TH2*) pfit->FindObject(pfit->GetName());
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
    if (doEps) pfit->Print(prefix+alias+".eps");
    c1->SetRightMargin(orm);
}

/** Reset line styles and colors, which get messed up by tdrStyle */
void prettyLine(TCanvas *canv, int pad, const char *cname, int color) {
    TGraph *c = (TGraph *) canv->GetPad(pad)->FindObject(cname);
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
            std::cerr << "Trying " << baff << std::endl;
            TDirectory *subdir = (TDirectory *) getFromPrefix(dir, baff, /*prefixOnly=*/false);
            if (subdir == 0) {
                //std::cerr << "Didn't find '" << baff << "*' in " << dir->GetName() << std::endl;
                continue;
            }
            subdir->ls();
            TCanvas *fitc = (TCanvas *) subdir->Get("fit_canvas");
            if (fitc == 0) {
                std::cerr << "Didn't find " << TString(baff) << "/fit_canvas in " << dir->GetName() << std::endl;
                continue;
            }
            fitc->Draw(); 
            prettyLines(fitc);
            fitc->Print(prefix+TString("canvases/")+buff+"_fit.png");
            fitc->Close();
        }
    }
}

enum GraphMergeAlgo { Combine, LowestUnc, Last };
TGraphAsymmErrors *mergeGraphs(TGraphAsymmErrors *g1, TGraphAsymmErrors *g2, TGraphAsymmErrors *g3=0, TGraphAsymmErrors *g4=0, GraphMergeAlgo algo=LowestUnc, double yMinCut = 0.001) {
    TGraphAsymmErrors *graphs[4];
    graphs[0] = g1;
    graphs[1] = g2;
    graphs[2] = g3;
    graphs[3] = g4;
    int ng = 0; int n = 0;
    TGraphAsymmErrors *gFirst = 0;
    for (int i = 0; i < 4; ++i) { 
        if (graphs[i] != 0 && graphs[i]->GetN() > 0) {
            ng++; n = TMath::Max(n, graphs[i]->GetN());
            if (gFirst == 0) gFirst = graphs[i];
        }
    }
    if (gFirst == 0) return 0;
    TGraphAsymmErrors *ret = new TGraphAsymmErrors(n);
    for (int i = 0; i < n; ++i) {
        double syw2 = 0, sw2 = 0, slo = 0, shi = 0;
        double xbin = gFirst->GetX()[i], sxw2 = 0;
        double xmax = xbin + gFirst->GetErrorXhigh(i), xmin = xbin - gFirst->GetErrorXlow(i);
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
            if (algo == Combine && ytry/ntry > yMinCut && y < yMinCut) continue;
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
        ret->SetPointError(i, sxw2/sw2 - xmin, xmax - sxw2/sw2, sqrt(slo)/sw2, sqrt(shi)/sw2);
    }
    return ret;
}    



/** Grab a TObject from a directory
 *  even knowing only the beginning 
 *  of it's name.
 */
TObject *getFromPrefix(TDirectory *dir, TString prefixName) {
    TIter next(dir->GetListOfKeys());
    TKey *key;
    while ((key = (TKey *) next())) {
        if (strstr(key->GetName(), prefixName.Data())) {
            return dir->Get(key->GetName());
        }
    }
    return 0;
}
void plotUtil() { }

/* Other useful plotting macros */

TString preliminary = ""; //"CMS Preliminary"
TString retitle = "";
TString datalbl = "Data", reflbl = "Sim.";

void cmsprelim() {
    TPaveText *cmsprel = new TPaveText(.70,.16,.94,.21,"NDC");
    cmsprel->SetTextSize(0.05);
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
    TLegend *leg = new TLegend(.78,.15 + legend_y_offset,.92,.27 + legend_y_offset);
    leg->AddEntry(g1, lab1, "LP");
    leg->AddEntry(g2, lab2, "LP");
    leg->SetTextSize(0.05);
    leg->SetTextFont(42);
    leg->SetShadowColor(0);
    leg->SetFillColor(0);
    if (preliminary != "") cmsprelim();
    leg->Draw();
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
void reTitleY(TCanvas *pl, TString ytitle) {
    TH1 *first = (TH1*) pl->GetListOfPrimitives()->At(0);
    TH1 *last = (TH1*) pl->GetListOfPrimitives()->At(pl->GetListOfPrimitives()->GetSize()-1);
    if (first) reTitleTAxis(first->GetYaxis(), ytitle);
    if (last)  reTitleTAxis(last->GetYaxis(), ytitle);
}   
void reTitleTAxis(TAxis *ax, TString ytitle, double yoffset=1.0) {
   ax->SetTitle(ytitle); 
   ax->SetTitleOffset(yoffset); 
   ax->SetDecimals(true);
}


void fixupErrors(TGraphAsymmErrors *gr) {
    for (size_t i = 0; i < gr->GetN(); ++i) {
        if (gr->GetErrorYhigh(i) == 0 && gr->GetErrorYlow(i) == 0) continue;
        if (gr->GetErrorYhigh(i) == 0 && gr->GetY()[i] != 1.0) {
            std::cerr << "Fixup error high for " << gr->GetName() << ", efficiency = " << gr->GetY()[i] << std::endl;
            gr->SetPointEYhigh(i, 1.0-gr->GetY()[i]);
        }
        if (gr->GetErrorYlow(i) == 0 && gr->GetY()[i] != 0.0) {
            std::cerr << "Fixup error low for " << gr->GetName() << ", efficiency = " << gr->GetY()[i] << std::endl;
            gr->SetPointEYlow(i, gr->GetY()[i]);
        }
    }
}

void doRatio(RooHist *hfit, RooHist *href, TString alias) {
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
    ratio.GetYaxis()->SetRangeUser(1-1.5*max,1+1.2*max);
    ratio.GetXaxis()->SetRangeUser(ratio.GetX()[0]-ratio.GetErrorXlow(0), ratio.GetX()[ratio.GetN()-1]+ratio.GetErrorXhigh(ratio.GetN()-1));
    if (datalbl) reTitleTAxis(ratio.GetYaxis(), datalbl+"/"+reflbl+" ratio");
    if (preliminary != "") cmsprelim();
    gPad->Print(prefix+alias+"_ratio.png");
}

void doDiff(RooHist *hfit, RooHist *href, TString alias) {
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
    diff.GetYaxis()->SetRangeUser(-1.5*max,1.2*max);
    if (datalbl) reTitleTAxis(diff.GetYaxis(), datalbl+" - "+reflbl+" difference");
    if (preliminary != "") cmsprelim();
    gPad->Print(prefix+alias+"_diff.png");

}
/** Plot FIT from file 1 plus FIT from file 2 */
void refstack(TDirectory *fit, TDirectory *ref, TString alias, TString fitname) {
    TCanvas *pref = getFromPrefix(ref->GetDirectory("fit_eff_plots"), fitname);
    if (pref == 0) {
        std::cerr << "NOT FOUND: " << "fit_eff_plots/"+fitname << " in " << ref->GetName() << std::endl;
        return;
    }
    RooHist *href = (RooHist *) pref->FindObject("hxy_fit_eff");
    fixupErrors(href);
    href->SetLineWidth(2);
    href->SetLineColor(kRed);
    href->SetMarkerColor(kRed);
    href->SetMarkerStyle(25);
    href->SetMarkerSize(2.0);

    TCanvas *pfit = getFromPrefix(fit->GetDirectory("fit_eff_plots"), fitname);
    if (pfit == 0) {
        std::cerr << "NOT FOUND: " << "fit_eff_plots/"+fitname << " in " << fit->GetName() << std::endl;
        return;
    }
    RooHist *hfit = (RooHist *) pfit->FindObject("hxy_fit_eff");
    fixupErrors(hfit);
    hfit->SetLineWidth(2);
    hfit->SetLineColor(kBlack);
    hfit->SetMarkerColor(kBlack);
    hfit->SetMarkerStyle(20);
    hfit->SetMarkerSize(1.6);

    if (retitle != "") reTitleY(pref, retitle);

    pref->Draw();
    hfit->Draw("P SAME");
    if (datalbl) doLegend(hfit,href,datalbl,reflbl);
    gPad->Print(prefix+alias+".png");

    doRatio(hfit,href,alias); 
    doDiff(hfit,href,alias); 
}

/** Plot FIT from file 1 plus CNT from file 2 */
void mcstack(TDirectory *fit, TDirectory *ref, TString alias, TString name) {
    TCanvas *pref = getFromPrefix(ref->GetDirectory("cnt_eff_plots"), name);
    RooHist *href = (RooHist *) pref->FindObject("hxy_cnt_eff");
    fixupErrors(href);
    href->SetLineWidth(2);
    href->SetLineColor(209);
    href->SetMarkerColor(209);
    href->SetMarkerStyle(25);
    href->SetMarkerSize(2.0);

    TCanvas *pfit = getFromPrefix(fit->GetDirectory("fit_eff_plots"), name);
    RooHist *hfit = (RooHist *) pfit->FindObject("hxy_fit_eff");
    fixupErrors(hfit);
    hfit->SetLineWidth(2);
    hfit->SetLineColor(206);
    hfit->SetMarkerColor(206);
    hfit->SetMarkerStyle(20);
    hfit->SetMarkerSize(1.6);

    if (retitle != "") reTitleY(pref, retitle);
    pref->Draw();
    hfit->Draw("P SAME");

    doLegend(hfit,href,datalbl,reflbl);
    gPad->Print(prefix+alias+".png");

    doRatio(hfit,href,alias); 
    doDiff(hfit,href,alias); 
}


/** Plot just one set */
void single( TDirectory *fit, TString alias, TString fitname) {
    TCanvas *pfit = getFromPrefix(fit->GetDirectory("fit_eff_plots"), fitname);
    if (pfit == 0) {
        std::cerr << "NOT FOUND: " << "fit_eff_plots/"+fitname << " in " << fit->GetName() << std::endl;
        return;
    }
    if (retitle != "") reTitleY(pfit, retitle);

    RooHist *hfit = (RooHist *) pfit->FindObject("hxy_fit_eff");
    fixupErrors(hfit);
    hfit->SetLineWidth(2);
    hfit->SetLineColor(kBlue);
    hfit->SetMarkerColor(kBlue);
    hfit->SetMarkerStyle(20);
    hfit->SetMarkerSize(1.6);
    if (preliminary != "") cmsprelim();
    pfit->Print(prefix+alias+".png"); 
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
            TDirectory *subdir = (TDirectory *) getFromPrefix(dir, baff);
            if (subdir == 0) {
                std::cerr << "Didn't find '" << baff << "*' in " << dir->GetName() << std::endl;
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

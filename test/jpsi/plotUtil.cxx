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

/** Plot FIT from file 1 plus FIT from file 2 */
void refstack(TDirectory *fit, TDirectory *ref, TString alias, TString fitname) {
    RooPlot *pref = getFromPrefix(ref->GetDirectory("fit_eff_plots"), fitname);
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

    RooPlot *pfit = getFromPrefix(fit->GetDirectory("fit_eff_plots"), fitname);
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

    pref->Draw();
    hfit->Draw("P SAME");
    gPad->Print(prefix+alias+".png");

    TGraphAsymmErrors diff(hfit->GetN()), zero(hfit->GetN());
    double max = 0;
    for (size_t i = 0, n = hfit->GetN(); i < n; ++i) {
        max = TMath::Max(max, fabs(hfit->GetY()[i] - href->GetY()[i]) + fabs(hfit->GetErrorYhigh(i)) + fabs(hfit->GetErrorYlow(i)));
        max = TMath::Max(max, fabs(href->GetErrorYlow(i)) + fabs(href->GetErrorYhigh(i)));
        diff.SetPoint(i, hfit->GetX()[i], hfit->GetY()[i] - href->GetY()[i]);
        diff.SetPointError(i, hfit->GetErrorXlow(i), hfit->GetErrorXhigh(i), 
                              hfit->GetErrorYlow(i), hfit->GetErrorYhigh(i));
        zero.SetPoint(i, href->GetX()[i], 0);
        zero.SetPointError(i, href->GetErrorXlow(i), href->GetErrorXhigh(i), 
                              href->GetErrorYlow(i), href->GetErrorYhigh(i));
    }
    zero.SetLineWidth(2);
    zero.SetLineColor(kRed);
    zero.SetMarkerColor(kRed);
    zero.SetMarkerStyle(25);
    zero.SetMarkerSize(2.0);
    diff.SetLineWidth(2);
    diff.SetLineColor(kBlack);
    diff.SetMarkerColor(kBlack);
    diff.SetMarkerStyle(20);
    diff.SetMarkerSize(1.6);

    //diff.Draw("AP");
    zero.Draw("AP"); // SAME");
    diff.Draw("P SAME");
    zero.GetYaxis()->SetRangeUser(-1.2*max,1.2*max);
    gPad->Print(prefix+alias+"_diff.png");
}

/** Plot FIT from file 1 plus CNT from file 2 */
void mcstack(TDirectory *fit, TDirectory *ref, TString alias, TString name) {
    RooPlot *pref = getFromPrefix(ref->GetDirectory("cnt_eff_plots"), name);
    RooHist *href = (RooHist *) pref->FindObject("hxy_cnt_eff");
    href->SetLineWidth(2);
    href->SetLineColor(209);
    href->SetMarkerColor(209);
    href->SetMarkerStyle(25);
    href->SetMarkerSize(2.0);

    RooPlot *pfit = getFromPrefix(fit->GetDirectory("fit_eff_plots"), name);
    RooHist *hfit = (RooHist *) pfit->FindObject("hxy_fit_eff");
    hfit->SetLineWidth(2);
    hfit->SetLineColor(206);
    hfit->SetMarkerColor(206);
    hfit->SetMarkerStyle(20);
    hfit->SetMarkerSize(1.6);

    pref->Draw();
    hfit->Draw("P SAME");
    gPad->Print(prefix+alias+".png");
}


/** Plot just one set */
void single( TDirectory *fit, TString alias, TString fitname) {
    RooPlot *pfit = getFromPrefix(fit->GetDirectory("fit_eff_plots"), fitname);
    if (pfit == 0) {
        std::cerr << "NOT FOUND: " << "fit_eff_plots/"+fitname << " in " << fit->GetName() << std::endl;
        return;
    }
    RooHist *hfit = (RooHist *) pfit->FindObject("hxy_fit_eff");
    hfit->SetLineWidth(2);
    hfit->SetLineColor(kBlue);
    hfit->SetMarkerColor(kBlue);
    hfit->SetMarkerStyle(20);
    hfit->SetMarkerSize(1.6);
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

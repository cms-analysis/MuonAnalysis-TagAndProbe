#include <TCanvas.h>
#include <TPad.h>
#include "plotUtil.cxx"
TString prefix = "plots_dev/muonid_collage/";
TString basedir  = "tpTree";

TFile *jpsiFile, *zmmFile;

TCanvas *c1 = 0;
TString year = "2012", energy = "8 TeV";
void collageMuonID_ZZ4L() {

    prefix = "plots/muonid_5fb/collage/";
    gSystem->mkdir(prefix,true);

    gROOT->ProcessLine(".x /afs/cern.ch/user/g/gpetrucc/cpp/tdrstyle.cc");
    gStyle->SetOptStat(0);
    c1 = new TCanvas("c1","c1");

    if (gROOT->GetListOfFiles()->GetEntries() == 2) {
        jpsiFile = (TFile*) gROOT->GetListOfFiles()->At(0);
        zmmFile  = (TFile*) gROOT->GetListOfFiles()->At(1);
        std::cout << "Will use " << jpsiFile->GetName() << " for J/Psi, " << zmmFile->GetName() << " for Z" << std::endl;
    } else {
        std::cout << "This macro requires two files (JPsi, Z)" << std::endl;
        return;
    }

    doFillMC = true;
    doRatioPlot = true;
    doDiffPlot = false;
    doPdf = true;
    doSquare = true; yMin = 0; yMax = 1.1;
    if (TString(gROOT->GetListOfFiles()->At(0)->GetName()).Contains("2011")) {
        year = "2011"; energy = "7 TeV";
        prefix += "2011_";
    } else {
        prefix += "2012_";
    }
    datalbl = "Data, "+year;
    reflbl  = "Simulation";
    //cmsprel_xoffs = 0.2;
    preliminary = "CMS Preliminary,  #sqrt{s} = "+energy;
    doLogX = true;

    yMinD = -0.179; yMaxD = 0.179;
    yMinR =  0.750; yMaxR = 1.25;

    fOut = TFile::Open(prefix+"fitComb.root", "UPDATE");

    plotMuonIDData();
}

void fixmeFixmeFixme(TGraphAsymmErrors *fit) {
    sanityCheck(fit);
    for (int i = 0; i < fit->GetN(); ++i) {
        double x = fit->GetX()[i], y = fit->GetY()[i], eyl = fit->GetErrorYlow(i), eyh = fit->GetErrorYhigh(i);
        if (doFillMC) {
            double errMin = (x < 20 ? 0.01 : 0.005);
            if (eyl <= errMin && y >   errMin) eyl = errMin;
            if (eyh <= errMin && y < 1-errMin) eyh = errMin;
        }
        fit->SetPointEYhigh(i, eyh);
        fit->SetPointEYlow(i, eyl);
    }
}


TGraphAsymmErrors * merge2(TGraphAsymmErrors *jpsiLo, TGraphAsymmErrors *jpsi, double ptMin, double ptLo, double ptZmm=999.) {
    TGraphAsymmErrors *ret = new TGraphAsymmErrors();
    if (jpsiLo != 0) {
        for (int i = 0, n = jpsiLo->GetN(); i < n; ++i) {
            double x = jpsiLo->GetX()[i];
            if (x < ptMin) continue;
            if (x <= ptLo) {
                int k = ret->GetN(); ret->Set(k+1);
                ret->SetPoint(k, jpsiLo->GetX()[i], jpsiLo->GetY()[i]);
                ret->SetPointError(k, jpsiLo->GetErrorXlow(i), jpsiLo->GetErrorXhigh(i), jpsiLo->GetErrorYlow(i), jpsiLo->GetErrorYhigh(i));
            } else break;
        }
    }
    if (jpsi != 0) {
        for (int i = 0, n = jpsi->GetN(); i < n; ++i) {
            double x = jpsi->GetX()[i];
            if ((x <= ptLo) && (jpsiLo != 0)) continue;
            if (x <= ptZmm) {
                int k = ret->GetN(); ret->Set(k+1);
                ret->SetPoint(k, jpsi->GetX()[i], jpsi->GetY()[i]);
                ret->SetPointError(k, jpsi->GetErrorXlow(i), jpsi->GetErrorXhigh(i), jpsi->GetErrorYlow(i), jpsi->GetErrorYhigh(i));
            } else break;
        }
    }
    if (ret->GetN() == 0) return NULL;
    fixmeFixmeFixme(ret);
    return ret;
}

TGraphAsymmErrors * merge3(TGraphAsymmErrors *jpsiLo, TGraphAsymmErrors *jpsi, TGraphAsymmErrors *zmm, double ptMin, double ptLo, double ptJPsi, double ptZmm=999.) {
    TGraphAsymmErrors *ret = new TGraphAsymmErrors();
    if (jpsiLo != 0) {
        for (int i = 0, n = jpsiLo->GetN(); i < n; ++i) {
            double x = jpsiLo->GetX()[i];
            if (x < ptMin) continue;
            if (x <= ptLo) {
                int k = ret->GetN(); ret->Set(k+1);
                ret->SetPoint(k, jpsiLo->GetX()[i], jpsiLo->GetY()[i]);
                ret->SetPointError(k, jpsiLo->GetErrorXlow(i), jpsiLo->GetErrorXhigh(i), jpsiLo->GetErrorYlow(i), jpsiLo->GetErrorYhigh(i));
            } else break;
        }
    }
    if (jpsi != 0) {
        for (int i = 0, n = jpsi->GetN(); i < n; ++i) {
            double x = jpsi->GetX()[i];
            if ((x <= ptLo) && (jpsiLo != 0)) continue;
            if (x <= ptJPsi || zmm == 0) {
                int k = ret->GetN(); ret->Set(k+1);
                ret->SetPoint(k, jpsi->GetX()[i], jpsi->GetY()[i]);
                ret->SetPointError(k, jpsi->GetErrorXlow(i), jpsi->GetErrorXhigh(i), jpsi->GetErrorYlow(i), jpsi->GetErrorYhigh(i));
            } else break;
        }
    }
    if (zmm != 0) {
        for (int i = 0, n = zmm->GetN(); i < n; ++i) {
            double x = zmm->GetX()[i];
            if (x <= ptJPsi) continue;
            if (x <= ptZmm) {
                int k = ret->GetN(); ret->Set(k+1);
                ret->SetPoint(k, zmm->GetX()[i], zmm->GetY()[i]);
                ret->SetPointError(k, zmm->GetErrorXlow(i), zmm->GetErrorXhigh(i), zmm->GetErrorYlow(i), zmm->GetErrorYhigh(i));
            } else break;
        }
    }
    if (ret->GetN() == 0) return NULL;
    fixmeFixmeFixme(ret);
    return ret;
}

TF1* paramTurnOn(TGraphAsymmErrors *graph, double xmin, double xmax, double xZ=999, bool twoTurnOn=false) {
    // TurnOn:  A + B*Fermi(x-xUp)) 
    // TurnOff: C * exp((x-x0)/T)
    TF1 *func = new TF1(Form("%s_func",graph->GetName()),"(x <= [0]) * ([1]/(1 + exp(-(x-[2])/[3]))+[4])*(1-0.1*pow(max(0,(20-x)/20),[5]))+ (x > [0]) * [6] * exp([7]*(x-[0])/[0])", xmin, xmax);
    func->SetParNames("switch", "plateuOn", "threshold", "resol", "baseline", "fix", "plateauOff", "decay");

    func->FixParameter(0, xZ);

    func->SetLineWidth(2);
    func->SetLineColor(kRed);

    double ymax = 0, ymaxZ = 0;
    int N = graph->GetN(); double *x = graph->GetX(), *y = graph->GetY();
    for (int i = 0; i < N; ++i) { 
        if (y[i] > ymax) ymax = y[i]; 
    }
    func->SetParameter(1, ymax);
    func->SetParLimits(1, ymax*0.8, ymax*1.2);
    func->SetParLimits(4, -0.2, 0.2);
    
    // wild guess for now 
    double xlo = 0, xhi = 0;
    for (int i = 0; i < N; ++i) {
        if (y[i] < 0.3*ymax) { xlo = x[i]; } else { break; } 
    }
    for (int i = N-1; i > 0; --i) {
        if (y[i] > 0.9*ymax) { xhi = x[i]; } else { break; } 
    }
    func->SetParameter(2, 0.5*(xhi+xlo));
    func->SetParLimits(2, xlo-1, xhi+1);
    func->SetParameter(3, 0.5*(xhi-xlo)+1);
    func->SetParLimits(3, 0.1, 5);

    func->SetParameter(5, 5);
    func->SetParLimits(5, 3., 7.);

    if (xZ < 100) {
        func->SetParameter(6, ymax);
        func->SetParLimits(6, ymax*0.8, ymax*1.2);
        func->SetParameter(7, 0.05);
        func->SetParLimits(7, -0.2, 0.2);
    } else {
        func->FixParameter(6, 0);
        func->FixParameter(7, 0);
    }
 
    graph->Fit(func,"B N W EX0","",xmin,xmax);
    
    return func;
}

TF1* paramTurnOnHLT(TGraphAsymmErrors *graph, double xmin, double xmax, double xZ=999, bool twoTurnOn=false) {
    // TurnOn:  B*Fermi(x-xUp)) 
    // TurnOff: B*exp((x-x0)/T)
    TF1 *func = new TF1(Form("%s_func",graph->GetName()),"(x <= [0]) * ([1]/(1 + exp(-(x-[2])/[3]))) +(x > [0]) * [1] * exp([4]*(x-[0])/[0])", xmin, xmax);
    func->SetParNames("switch", "plateu", "threshold", "resol", "decay");

    func->FixParameter(0, xZ);

    func->SetLineWidth(2);
    func->SetLineColor(kRed);

    double ymax = 0, ymaxZ = 0;
    int N = graph->GetN(); double *x = graph->GetX(), *y = graph->GetY();
    for (int i = 0; i < N; ++i) { 
        if (y[i] > ymax) ymax = y[i]; 
    }
    func->SetParameter(1, ymax);
    func->SetParLimits(1, ymax*0.8, ymax*1.2);
    
    // wild guess for now 
    double xlo = 0, xhi = 0;
    for (int i = 0; i < N; ++i) {
        if (y[i] < 0.3*ymax) { xlo = x[i]; } else { break; } 
    }
    for (int i = N-1; i > 0; --i) {
        if (y[i] > 0.9*ymax) { xhi = x[i]; } else { break; } 
    }
    func->SetParameter(2, 0.5*(xhi+xlo));
    func->SetParLimits(2, xlo-1, xhi+1);
    func->SetParameter(3, 0.5*(xhi-xlo)+0.5);
    func->SetParLimits(3, 0.01, 1);

    if (xZ < 100) {
        func->SetParameter(4, 0.05);
        func->SetParLimits(4, -0.2, 0.2);
    } else {
        func->FixParameter(4, 0);
    }
 
    graph->Fit(func,"B N W EX0","",xmin,xmax);
    
    return func;
}



void plotMuonIDData() {
    retitle = "Efficiency";

    const int nids  =  18;
    char *ids[nids]    = {  "PF",  "Glb",   "Tight2012noIP",  "TMA",  "ZZLoose" ,   
                            "SIP4_from_PF", "PFIso40_from_PF_and_SIP4",
                            "PFIso20DB_from_PF_and_SIP4", "PFIso12DB_from_PF_and_SIP4",
                            "PFIso40_from_Tight2012withIP", "PFIso20DB_from_Tight2012withIP", "PFIso12DB_from_Tight2012withIP",
                            "HLT_Mu17_from_PF_and_SIP4_and_PFIso40", "HLT_Mu8_from_PF_and_SIP4_and_PFIso40", "HLT_OrMu17_from_PF_and_SIP4_and_PFIso40", "HLT_OrMu8_from_PF_and_SIP4_and_PFIso40", "HLT_TkMu17_from_PF_and_SIP4_and_PFIso40", "HLT_TkMu8_from_PF_and_SIP4_and_PFIso40"  };
    char *titles[nids] = { "PF ID", "Glb ID",   "Tight ID (no ip)",   "TM Arb ID", "ZZ Loose ID",  
                           "SIP < 4",      "PF Iso < 0.4" ,          
                           "PF Rel Iso < 0.2 (#Delta#beta)", "PF Rel Iso < 0.12 (#Delta#beta)",
                           "PF Rel Iso < 0.4", "PF Rel Iso < 0.2 (#Delta#beta)", ,"PF Rel Iso < 0.12 (#Delta#beta)",
                           "Mu17",     "Mu8", "Mu17 OR",     "Mu8 OR", "Mu17 Tk", "Mu8 Tk" };

    const int neta  = 2;
    char *eta[neta]    = { "barrel", "endcaps" };
    char *teta[neta]   = { "        |#eta| < 1.2", "  1.2 < |#eta| < 2.4" };

    for (size_t i = 0; i < nids; ++i) {
        //if (i != 0) continue;
        for (size_t j = 0; j < neta; ++j) {
            for (size_t binning = 1; binning <= 2; ++binning) {
                //if (binning != 2) continue;
                yMin = 0.0; yMax = 1.1;
                TString idname(ids[i]);
                retitle = TString(titles[i])+" efficiency";
                TString etaname(eta[j]);

                TString jpsiBinName = ""; if (binning == 2) jpsiBinName = "2";
                TString zmmBinName = "";  if (binning == 2) zmmBinName = "2";

                TGraphAsymmErrors *fitJPsi, *fitZmm, *fitComb;
                TGraphAsymmErrors *refJPsi, *refZmm, *refComb;
                fitZmm = (TGraphAsymmErrors *) zmmFile->Get("fit_"+idname+"_pt_"+etaname+zmmBinName);
                refZmm = (TGraphAsymmErrors *) zmmFile->Get("ref_"+idname+"_pt_"+etaname+zmmBinName);
                if (fitZmm == 0 || refZmm == 0) continue;
                double ptMin = (binning == 1 ? 2.99 : 4.99), zPt = 15., ptMax = 999.;
                if (idname.Contains("Mu8")) ptMin = 5;
                if (idname.Contains("Mu17")) ptMin = 10;
                //if (idname.Contains("Mu17")) { ptMin = 15; ptMax = 40; }
                if (year == "2011") {
                    TGraphAsymmErrors *fitJPsi2 = (TGraphAsymmErrors *) jpsiFile->Get("fit_"+idname+"_pt_Mu5_Track2_"+etaname+jpsiBinName);
                    TGraphAsymmErrors *refJPsi2 = (TGraphAsymmErrors *) jpsiFile->Get("ref_"+idname+"_pt_Mu5_Track2_"+etaname+jpsiBinName);
                    TGraphAsymmErrors *fitJPsi7 = (TGraphAsymmErrors *) jpsiFile->Get("fit_"+idname+"_pt_Mu7_Track7_"+etaname+jpsiBinName);
                    TGraphAsymmErrors *refJPsi7 = (TGraphAsymmErrors *) jpsiFile->Get("ref_"+idname+"_pt_Mu7_Track7_"+etaname+jpsiBinName);
                    fitJPsi = merge2(fitJPsi2, fitJPsi7, 2.0, 7.5);
                    refJPsi = merge2(refJPsi2, refJPsi7, 2.0, 7.5);
                } else {
                    TGraphAsymmErrors *fitJPsi2   = (TGraphAsymmErrors *) jpsiFile->Get("fit_"+idname+"_pt_Mu5_Track2_"  +etaname+jpsiBinName);
                    TGraphAsymmErrors *refJPsi2   = (TGraphAsymmErrors *) jpsiFile->Get("ref_"+idname+"_pt_Mu5_Track2_"  +etaname+jpsiBinName);
                    TGraphAsymmErrors *fitJPsi3p5 = (TGraphAsymmErrors *) jpsiFile->Get("fit_"+idname+"_pt_Mu5_Track3p5_"+etaname+jpsiBinName);
                    TGraphAsymmErrors *refJPsi3p5 = (TGraphAsymmErrors *) jpsiFile->Get("ref_"+idname+"_pt_Mu5_Track3p5_"+etaname+jpsiBinName);
                    //TGraphAsymmErrors *refJPsi3p5 = fitJPsi2; // FIXME wrong 2012 sample
                    TGraphAsymmErrors *fitJPsi7   = (TGraphAsymmErrors *) jpsiFile->Get("fit_"+idname+"_pt_Mu7_Track7_"  +etaname+jpsiBinName);
                    TGraphAsymmErrors *refJPsi7   = (TGraphAsymmErrors *) jpsiFile->Get("ref_"+idname+"_pt_Mu7_Track7_"  +etaname+jpsiBinName);
                    fitJPsi = merge3(fitJPsi2, fitJPsi3p5, fitJPsi7, 2.0, 3.75, 7.5);
                    refJPsi = merge3(refJPsi2, refJPsi3p5, refJPsi7, 2.0, 3.75, 7.5);
                }
                TGraphAsymmErrors *fit = merge2(fitJPsi, fitZmm, ptMin, zPt, ptMax);
                TGraphAsymmErrors *ref = merge2(refJPsi, refZmm, ptMin, zPt, ptMax);

                if (doFillMC) {
                    ref->SetLineWidth(1);
                    ref->SetLineColor(62);
                    ref->SetLineStyle(0);
                    ref->SetFillColor(65);
                    ref->SetMarkerColor(2);
                    ref->SetMarkerStyle(1);
                    ref->SetMarkerSize(0);
                } else {
                    ref->SetLineWidth(2);
                    ref->SetLineColor(kRed);
                    ref->SetMarkerColor(kRed);
                    ref->SetMarkerStyle(25);
                    ref->SetMarkerSize(2.0);
                }

                fit->SetLineWidth(2);
                fit->SetLineColor(kBlack);
                fit->SetMarkerColor(kBlack);
                fit->SetMarkerStyle(20);
                fit->SetMarkerSize(1.2);

                c1->cd(); c1->Clear(); c1->SetLogx(1); 
                if (doSquare) squareCanvas(c1);
                double xmax = TMath::Max(xmaxGraph(fit), xmaxGraph(ref));
                TH1F *frame = new TH1F("frame","frame",1,ptMin-0.01,xmax); 
                frame->GetYaxis()->SetRangeUser(yMin,yMax);
                frame->GetYaxis()->SetTitle(retitle);
                frame->GetXaxis()->SetTitle("muon p_{T}  (GeV/c)");
                frame->GetXaxis()->SetMoreLogLabels(1);
                frame->GetXaxis()->SetNoExponent(1);

                frame->Draw();
                TLine line; line.SetLineWidth(2.0); line.SetLineStyle(2); line.SetLineColor(202);
                if (fitJPsi != 0 && xmax > zPt) line.DrawLine(zPt, 0.5, zPt, yMax);

                ref->Draw(doFillMC ? "E2" : "P");
                fit->Draw("P");
                extraSpam = teta[j]; doLegend(fit, ref, datalbl, reflbl);

                TString binName = ""; if (binning == 2) binName = "2";
                c1->Print(prefix+idname+"_"+etaname+binName+".eps");
                gSystem->Exec("convert "+prefix+idname+"_"+etaname+binName+".eps "+prefix+idname+"_"+etaname+binName+".png");
                if (doPdf) c1->Print(prefix+idname+"_"+etaname+binName+".pdf");
                if (doEps) c1->Print(prefix+idname+"_"+etaname+binName+".eps");
                if (doTxt) printGraphs(fit,ref,idname+"_"+etaname+binName);

                if (fOut) fOut->WriteTObject(fit, "fit_"+idname+"_"+etaname, "Overwrite");
                if (fOut) fOut->WriteTObject(ref, "ref_"+idname+"_"+etaname, "Overwrite");

                if (idname.Contains("HLT")) {
                    TF1 *fitFunc = paramTurnOnHLT(fit, ptMin+0.5, xmax, 200); 
                    TF1 *refFunc = paramTurnOnHLT(ref, ptMin+0.5, xmax, 200); 
                    fitFunc->SetLineColor(100); fitFunc->SetLineWidth(2);
                    refFunc->SetLineColor(214); refFunc->SetLineWidth(2);

                    frame->Draw();
                    TLine line; line.SetLineWidth(2.0); line.SetLineStyle(2); line.SetLineColor(202);
                    if (idname.Contains("17")) line.DrawLine(17, 0.0, 17, yMax);
                    if (idname.Contains("8")) line.DrawLine(8, 0.0, 8, yMax);
                    ref->Draw(doFillMC ? "E2" : "P");
                    fit->Draw("P");
                    refFunc->Draw("L SAME");
                    fitFunc->Draw("L SAME");
                    doLegend(fit, ref, datalbl, reflbl);
                    c1->Print(prefix+idname+"_"+etaname+binName+"_func.eps");
                    gSystem->Exec("convert "+prefix+idname+"_"+etaname+binName+"_func.eps "+prefix+idname+"_"+etaname+binName+"_func.png");
                    if (doPdf) c1->Print(prefix+idname+"_"+etaname+binName+"_func.pdf");
                }

                TString databk = datalbl, refbk = reflbl;
                datalbl = "Data"; reflbl = "Sim.";
                if (doRatioPlot) doRatio(fit,ref,idname+"_"+etaname+binName,"muon p_{T}  (GeV/c)");
                if (doDiffPlot)   doDiff(fit,ref,idname+"_"+etaname+binName,"muon p_{T}  (GeV/c)");
                datalbl = databk; reflbl = refbk; 
                delete frame;

            }
        }
    }
}


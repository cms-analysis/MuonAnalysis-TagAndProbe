#include <TCanvas.h>
#include <TPad.h>
#include "plotUtil.cxx"

TFile *jpsiFile[2], *zmmFile[2];

void overallMuonID_ZZ4L() {
    jpsiFile[0] = TFile::Open("plots/muonid_v2/jpsi/2011_fitJPsi.root");
    jpsiFile[1] = TFile::Open("plots/muonid_v2/jpsi/2012_fitJPsi.root");
    zmmFile[0]  = TFile::Open("plots/muonid_v2/zmumu/2011_fitZ.root");
    zmmFile[1]  = TFile::Open("plots/muonid_v2/zmumu/2012_fitZ.root");

    fOut = TFile::Open("final.root", "RECREATE");
    
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

void makeSF(TDirectory *out, TString outname, TGraphAsymmErrors *data, TGraphAsymmErrors *ref) {
    TGraphAsymmErrors *ret = new TGraphAsymmErrors();
    for (int i = 0, k = 0; i < ref->GetN(); ++i) {
        int j = findBin(data,ref->GetX()[i]); if (j == -1) continue ;
        if (fabs(ref->GetY()[j]) < 0.05) continue; else ++k;
        double r   = data->GetY()[i]/ref->GetY()[j];
        double rup = (data->GetY()[i] == 0 ? data->GetErrorYhigh(i)/(ref->GetY()[j]) :
                                             r*TMath::Hypot(data->GetErrorYhigh(i)/data->GetY()[i], ref->GetErrorYlow(j)/ref->GetY()[j]));
        double rdn = (data->GetY()[i] == 0 ? 0 :
                                             r*TMath::Hypot(data->GetErrorYlow(i)/data->GetY()[i],  ref->GetErrorYhigh(j)/ref->GetY()[j]));
        ret->Set(k);
        ret->SetPoint(k-1, data->GetX()[i], r);
        ret->SetPointError(k-1, data->GetErrorXlow(i), data->GetErrorXhigh(i), rdn, rup);
    }
    
    ret->SetName(outname);
    out->WriteTObject(ret, outname);
    out->WriteTObject(data->Clone(outname+"_data"), outname+"_data");
    out->WriteTObject(ref->Clone(outname+"_ref"), outname+"_ref");
}
void makeSF(TDirectory *out, TString outname, TDirectory *in, TString inname) {
    TGraphAsymmErrors *data = (TGraphAsymmErrors *) in->Get("fit_"+inname); 
    if (data == 0) { std::cerr << "ERROR: missing fit_" << inname << " in " << in << "\n"; in->ls(); exit(1); }
    TGraphAsymmErrors *ref  = (TGraphAsymmErrors *) in->Get("ref_"+inname);
    if (ref == 0) { std::cerr << "ERROR: missing ref_" << inname << " in " << in << "\n"; in->ls(); exit(1); }
    makeSF(out, outname, data, ref);
}

void makeTH2F(TDirectory *out, TString outname, TString ptplot, TString etaplot) {
    TGraphAsymmErrors *effptb = (TGraphAsymmErrors *) out->Get(Form(ptplot.Data(), "barrel"));
    TGraphAsymmErrors *effpte = (TGraphAsymmErrors *) out->Get(Form(ptplot.Data(), "endcaps"));
    TGraphAsymmErrors *effeta = (TGraphAsymmErrors *) out->Get(etaplot);
    const int neta = 16;
    double etabins[neta] = { -2.4, -2.1, -1.6, -1.2, -0.9, -0.6, -0.3, -0.2, 0.2, 0.3, 0.6, 0.9, 1.2, 1.6, 2.1, 2.4 };
    const int npt = 8;
    double ptbins[npt] = { 5, 7.5, 10, 15, 20, 30, 40, 100 };
    TH2F *ret = new TH2F(outname,outname, npt-1, ptbins, neta-1, etabins);
    for (int i = 1; i < npt; ++i) {
        for (int j = 1; j < neta; ++j) {
            double pt = 0.5*(ptbins[i] + ptbins[i-1]), eta = 0.5*(etabins[j] + etabins[j-1]);
            if (pt > 20) {
                int k = findBin(effeta, eta);
                if (k == -1) { std::cerr << "ERROR, couldn't find point at eta " << eta << " in " << etaplot << std::endl; continue; }
                ret->SetBinContent(i, j, effeta->GetY()[k]);
                ret->SetBinError(i, j,  TMath::Max(effeta->GetErrorYhigh(k), effeta->GetErrorYlow(k)));;
            } else {
                TGraphAsymmErrors *effpt = (fabs(eta) < 1.2 ? effptb : effpte);
                int k = findBin(effpt, pt);
                if (k == -1) { std::cerr << "ERROR, couldn't find point at pt " << pt << ", eta " << eta << " in " << ptplot << std::endl; continue; }
                ret->SetBinContent(i, j, effpt->GetY()[k]);
                ret->SetBinError(i, j,  TMath::Max(effpt->GetErrorYhigh(k), effpt->GetErrorYlow(k)));;
            }
        }
    }
    out->WriteTObject(ret,outname);
}

void plotMuonIDData() {
    // Muon ID, high pT
    makeSF(fOut, "ID_eta_2011A", zmmFile[0], "PF_from_2011A_eta");
    makeSF(fOut, "ID_eta_2011B", zmmFile[0], "PF_from_2011B_eta");
    makeSF(fOut, "ID_eta_2012",  zmmFile[1], "PF_eta");

    // Muon SIP, high pT
    makeSF(fOut, "IP_eta_2011", zmmFile[0], "SIP4_from_PF_eta");
    makeSF(fOut, "IP_eta_2012", zmmFile[1], "SIP4_from_PF_eta");

    // Muon ISO, high pT (useless)
    makeSF(fOut, "ISO_eta_2011", zmmFile[0], "PFIso40_from_PF_and_SIP4_eta");
    makeSF(fOut, "ISO_eta_2012", zmmFile[1], "PFIso40_from_PF_and_SIP4_eta");

    // Muon Trigger, high pT
    makeSF(fOut, "Mu17_eta_2011A", zmmFile[0], "HLT_Mu17_from_PF_and_SIP4_and_PFIso40_and_2011A_eta");
    makeSF(fOut, "Mu17_eta_2011B", zmmFile[0], "HLT_Mu17_from_PF_and_SIP4_and_PFIso40_and_2011B_eta");
    makeSF(fOut, "Mu17_eta_2012",  zmmFile[1], "HLT_OrMu17_from_PF_and_SIP4_and_PFIso40_eta");
    makeSF(fOut, "Mu8_eta_2011A", zmmFile[0], "HLT_Mu8_from_PF_and_SIP4_and_PFIso40_and_2011A_eta");
    makeSF(fOut, "Mu8_eta_2011B", zmmFile[0], "HLT_Mu8_from_PF_and_SIP4_and_PFIso40_and_2011B_eta");
    makeSF(fOut, "Mu8_eta_2012",  zmmFile[1], "HLT_OrMu8_from_PF_and_SIP4_and_PFIso40_eta");

    // Muon ID, low pT
    for (int y = 0; y <= 1; ++y) {
        for (int i = 0; i < 2; ++i) {
            TString etaname = (i == 0 ? "barrel" : "endcaps");
            TGraphAsymmErrors *fitJPsi1, *fitJPsi2, *fitZmm, *fitComb;
            TGraphAsymmErrors *refJPsi1, *refJPsi2, *refZmm, *refComb;
            fitZmm = (TGraphAsymmErrors *) zmmFile[y]->Get("fit_PF_pt_"+etaname);
            refZmm = (TGraphAsymmErrors *) zmmFile[y]->Get("ref_PF_pt_"+etaname);
            fitJPsi1 = (TGraphAsymmErrors *) jpsiFile[y]->Get("fit_PF"+TString(y ? "_pt_Mu5_Track3p5_" : "_pt_Mu5_Track2_")+etaname+"2");
            fitJPsi2 = (TGraphAsymmErrors *) jpsiFile[y]->Get("fit_PF_pt_Mu7_Track7_"+etaname+"2");
            refJPsi1 = (TGraphAsymmErrors *) jpsiFile[y]->Get("ref_PF"+TString(y ? "_pt_Mu5_Track3p5_" : "_pt_Mu5_Track2_")+etaname+"2");
            refJPsi2 = (TGraphAsymmErrors *) jpsiFile[y]->Get("ref_PF_pt_Mu7_Track7_"+etaname+"2");
            fitComb = merge3(fitJPsi1, fitJPsi2, fitZmm, 5, 7.5, 15, 20);
            refComb = merge3(refJPsi1, refJPsi2, refZmm, 5, 7.5, 15, 20);
            makeSF(fOut, Form("ID_pt_%s_%d", etaname.Data(), 2011+y), fitComb, refComb);
        }
    }
    makeTH2F(fOut, "TH2D_ID_2011A", "ID_pt_%s_2011", "ID_eta_2011A");
    makeTH2F(fOut, "TH2D_ID_2011B", "ID_pt_%s_2011", "ID_eta_2011B");
    makeTH2F(fOut, "TH2D_ID_2012",  "ID_pt_%s_2012", "ID_eta_2012");

    // Muon IP, low pT
    for (int y = 0; y <= 1; ++y) {
        for (int i = 0; i < 2; ++i) {
            TString etaname = (i == 0 ? "barrel" : "endcaps");
            TGraphAsymmErrors *fitZmm, *refZmm, *fitCrop, *refCrop;
            fitZmm = (TGraphAsymmErrors *) zmmFile[y]->Get("fit_SIP4_from_PF_pt_"+etaname+"2");
            refZmm = (TGraphAsymmErrors *) zmmFile[y]->Get("ref_SIP4_from_PF_pt_"+etaname+"2");
            if (fitZmm == 0) { std::cerr << "ERROR: missing fit_SIP4_from_PF_pt_"+etaname+"2 in " << in << "\n"; zmmFile[y]->ls(); exit(1); }
            if (refZmm == 0) { std::cerr << "ERROR: missing ref_SIP4_from_PF_pt_"+etaname+"2 in " << in << "\n"; zmmFile[y]->ls(); exit(1); }
            fitCrop = merge2(fitZmm, fitZmm, 5., 5., 20.);
            refCrop = merge2(refZmm, refZmm, 5., 5., 20.);
            makeSF(fOut, Form("IP_pt_%s_%d", etaname.Data(), 2011+y), fitCrop, refCrop);
        }
    }
    makeTH2F(fOut, "TH2D_IP_2011", "IP_pt_%s_2011", "IP_eta_2011");
    makeTH2F(fOut, "TH2D_IP_2012", "IP_pt_%s_2012", "IP_eta_2012");

    // Muon ISO, low pT
    for (int y = 0; y <= 1; ++y) {
        for (int i = 0; i < 2; ++i) {
            TString etaname = (i == 0 ? "barrel" : "endcaps");
            TGraphAsymmErrors *fitZmm, *refZmm, *fitCrop, *refCrop;
            fitZmm = (TGraphAsymmErrors *) zmmFile[y]->Get("fit_PFIso40_from_PF_and_SIP4_pt_"+etaname+"2");
            refZmm = (TGraphAsymmErrors *) zmmFile[y]->Get("ref_PFIso40_from_PF_and_SIP4_pt_"+etaname+"2");
            if (fitZmm == 0) { std::cerr << "ERROR: missing fit_PFIso40_from_PF_and_SIP4_pt_"+etaname+"2 in " << in << "\n"; zmmFile[y]->ls(); exit(1); }
            if (refZmm == 0) { std::cerr << "ERROR: missing ref_PFIso40_from_PF_and_SIP4_pt_"+etaname+"2 in " << in << "\n"; zmmFile[y]->ls(); exit(1); }
            fitCrop = merge2(fitZmm, fitZmm, 5., 5., 20.);
            refCrop = merge2(refZmm, refZmm, 5., 5., 20.);
            makeSF(fOut, Form("ISO_pt_%s_%d", etaname.Data(), 2011+y), fitCrop, refCrop);
        }
    }

    makeTH2F(fOut, "TH2D_ISO_2011", "ISO_pt_%s_2011", "ISO_eta_2011");
    makeTH2F(fOut, "TH2D_ISO_2012", "ISO_pt_%s_2012", "ISO_eta_2012");

    TH2F *hAll2011A = (TH2F*) fOut->Get("TH2D_ID_2011A")->Clone("TH2D_ALL_2011A");
    hAll2011A->Multiply((TH2F*) fOut->Get("TH2D_IP_2011"));
    hAll2011A->Multiply((TH2F*) fOut->Get("TH2D_ISO_2011"));
    fOut->WriteTObject(hAll2011A, "TH2D_ALL_2011A");
    TH2F *hAll2011B = (TH2F*) fOut->Get("TH2D_ID_2011B")->Clone("TH2D_ALL_2011B");
    hAll2011B->Multiply((TH2F*) fOut->Get("TH2D_IP_2011"));
    hAll2011B->Multiply((TH2F*) fOut->Get("TH2D_ISO_2011"));
    fOut->WriteTObject(hAll2011B, "TH2D_ALL_2011B");
    TH2F *hAll2012 = (TH2F*) fOut->Get("TH2D_ID_2012")->Clone("TH2D_ALL_2012");
    hAll2012->Multiply((TH2F*) fOut->Get("TH2D_IP_2012"));
    hAll2012->Multiply((TH2F*) fOut->Get("TH2D_ISO_2012"));
    fOut->WriteTObject(hAll2012, "TH2D_ALL_2012");

}


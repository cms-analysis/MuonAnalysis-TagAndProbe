#include <TCanvas.h>
#include <TPad.h>
//TCanvas *c1 = 0;

#include "plotUtil.cxx"

bool doCanvases = false;

TGraphAsymmErrors *merge(TGraphAsymmErrors **fits, int ntrig, int neta, int ieta) ;
void plotMuonIDData() ;

void plotMuonID_JPsi_ZZ4L(TString fOutName="") {
    basedir  = "tpTree";
    prefix = "plots/muonid_5fb/jpsi/";

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
    doPdf = true;
    doSquare = true; yMin = 0; yMax = 1.048;
    doFillMC = true;
    TString year = "2012", energy = "8 TeV";
    if (TString(gROOT->GetListOfFiles()->At(0)->GetName()).Contains("pt_abseta")) {
        fOutName = "fitJPsi.root";
    }
    if (TString(gROOT->GetListOfFiles()->At(0)->GetName()).Contains("2011")) {
        year = "2011"; energy = "7 TeV";
        if (TString(gROOT->GetListOfFiles()->At(0)->GetName()).Contains("2011A")) {
            prefix += "2011A_";
        } else if (TString(gROOT->GetListOfFiles()->At(0)->GetName()).Contains("2011B")) {
            prefix += "2011B_";
        } else {
            prefix += "2011_";
        }
    } else {
        prefix += "2012_";
    }
    datalbl = "Data, "+year;
    reflbl  = "Simulation";
    preliminary = "CMS Preliminary,  #sqrt{s} = "+energy;


    if (fOutName != "") fOut = TFile::Open(prefix+fOutName, "UPDATE");
    ((TFile*) gROOT->GetListOfFiles()->At(0))->cd();

    doCanvases = false; //if (scenario.Contains("data_vs_mc")) doCanvases = true;

    plotMuonIDData();
}

void plotMuonIDData() {
    retitle = "Efficiency";

    const int nids  = 6;
    const char *ids[nids]    = { "PF", "Tight2012noIP", "Tight2012withIP", "ZZLoose",  "Glb", "TMA"     };
    const char *titles[nids] = { "PF", "Tight (no IP)", "Tight (with IP)", "ZZ Loose", "Glb", "TM Arb"  };

   
    for (int i = 0; i < nids; ++i) {
        TString idname(ids[i]);
        retitle = TString(titles[i])+" muon efficiency";

        const int ntrig = 3;
        const char *trigname[ntrig] = { "Mu7_Track7", "Mu5_Track3p5", "Mu5_Track2" };

        for (int j = 0; j < ntrig; ++j) {
            TString trig = trigname[j];

            TDirectory *fit_pt_eta = gFile->GetDirectory(basedir+"/"+idname+"_pt_abseta_"+trig+"/");
            if (fit_pt_eta) {
                yMin = 0.0; yMax = 1.1;
                if (ref != 0) {
                    TDirectory *ref_pt_eta = ref->GetDirectory(basedir+"/"+idname+"_pt_abseta_"+trig+"/");
                    if (ref_pt_eta == 0 && trig == "Mu5_Track3p5") {
                        ref_pt_eta = ref->GetDirectory(basedir+"/"+idname+"_pt_abseta_"+"Mu5_Track2"+"/");
                    }
                    extraSpam = "        |#eta| < 1.2"; refstack(fit_pt_eta, ref_pt_eta, idname+"_pt_"+trig+"_barrel",  "pt_PLOT_abseta_bin0");
                    extraSpam = "  1.2 < |#eta| < 2.4"; refstack(fit_pt_eta, ref_pt_eta, idname+"_pt_"+trig+"_endcaps", "pt_PLOT_abseta_bin1");
                } else {
                    extraSpam = "        |#eta| < 1.2"; single(fit_pt_eta, idname+"_pt_"+trig+"_barrel",  "pt_PLOT_abseta_bin0");
                    extraSpam = "  1.2 < |#eta| < 2.4"; single(fit_pt_eta, idname+"_pt_"+trig+"_endcaps", "pt_PLOT_abseta_bin1");
                }
            }

            TDirectory *fit_pt_eta = gFile->GetDirectory(basedir+"/"+idname+"_pt_abseta2_"+trig+"/");
            if (fit_pt_eta) {
                yMin = 0.0; yMax = 1.1;
                if (ref != 0) {
                    TDirectory *ref_pt_eta = ref->GetDirectory(basedir+"/"+idname+"_pt_abseta2_"+trig+"/");
                    if (ref_pt_eta == 0 && trig == "Mu5_Track3p5") {
                        ref_pt_eta = ref->GetDirectory(basedir+"/"+idname+"_pt_abseta2_"+"Mu5_Track2"+"/");
                    }
                    extraSpam = "        |#eta| < 1.2"; refstack(fit_pt_eta, ref_pt_eta, idname+"_pt_"+trig+"_barrel2",  "pt_PLOT_abseta_bin0");
                    extraSpam = "  1.2 < |#eta| < 2.4"; refstack(fit_pt_eta, ref_pt_eta, idname+"_pt_"+trig+"_endcaps2", "pt_PLOT_abseta_bin1");
                } else {
                    extraSpam = "        |#eta| < 1.2"; single(fit_pt_eta, idname+"_pt_"+trig+"_barrel2",  "pt_PLOT_abseta_bin0");
                    extraSpam = "  1.2 < |#eta| < 2.4"; single(fit_pt_eta, idname+"_pt_"+trig+"_endcaps2", "pt_PLOT_abseta_bin1");
                }
            }


        }
    }

}


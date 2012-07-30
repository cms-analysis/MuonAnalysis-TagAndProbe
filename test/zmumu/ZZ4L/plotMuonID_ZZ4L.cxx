#include <TCanvas.h>
#include <TPad.h>
#include "plotUtil.cxx"

void plotMuonIDData() ;

void plotMuonID_ZZ4L(TString fOutName="") {

    prefix = "plots/muonid_5fb/zmumu/";
    gSystem->mkdir(prefix,true);

    gROOT->ProcessLine(".x tdrstyle.cc");
    gStyle->SetOptStat(0);
    c1 = new TCanvas("c1","c1");

    if (gROOT->GetListOfFiles()->GetEntries() == 2) {
        ref = (TFile *) gROOT->GetListOfFiles()->At(1);
        ((TFile*) gROOT->GetListOfFiles()->At(0))->cd();
    }

    doRatioPlot = true; yMaxR = 1.25; yMinR = 0.75;
    doDiffPlot = false;
    doPdf = true;
    doSanity = true;
    doSquare = true; yMin = 0.0; yMax = 1.048;
    doFillMC = true;
    TString year = "2012", energy = "8 TeV";
    fOutName = "fitZ.root";
    if (TString(gROOT->GetListOfFiles()->At(0)->GetName()).Contains("2011")) {
        year = "2011"; energy = "7 TeV";
        prefix += "2011_";
    } else {
        prefix += "2012_";
    }
    datalbl = "Data, "+year;
    reflbl  = "Simulation";
    preliminary = "CMS Preliminary,  #sqrt{s} = "+energy;

    if (fOutName != "") fOut = TFile::Open(prefix+fOutName, "UPDATE");
    ((TFile*) gROOT->GetListOfFiles()->At(0))->cd();

    plotMuonIDData();
}

void plotMuonIDData() {
    retitle = "Efficiency";

    const int nids  = 34;
    const char *ids[nids]    = { "PF", "SIP4_from_PF", "PFIso40_from_PF_and_SIP4" ,
                                 ,"PF_from_2011A", "PF_from_2011B"
                                 ,"Glb", "TMA", "ZZLoose", "Tight2012noIP", "Tight2012withIP"
                                 ,"PFIso20DB_from_PF_and_SIP4", "PFIso12DB_from_PF_and_SIP4"
                                 ,"PFIso40_from_Tight2012withIP", "PFIso20DB_from_Tight2012withIP", "PFIso12DB_from_Tight2012withIP"
                                 ,"HLT_Mu17_from_PF_and_SIP4_and_PFIso40", "HLT_Mu8_from_PF_and_SIP4_and_PFIso40"
                                 ,"DoubleMu17Mu8_Mu17_from_PF_and_SIP4_and_PFIso40", "DoubleMu17Mu8_Mu8_from_PF_and_SIP4_and_PFIso40"
                                 ,"DoubleMu17TkMu8_Mu17_from_PF_and_SIP4_and_PFIso40", "DoubleMu17TkMu8_TkMu8_from_PF_and_SIP4_and_PFIso40" 
                                 ,"DoubleMu17Mu8_Mu17_from_PF_and_SIP4_and_PFIso40_and_PostTS", "DoubleMu17Mu8_Mu8_from_PF_and_SIP4_and_PFIso40_and_PostTS"
                                 ,"DoubleMu17TkMu8_Mu17_from_PF_and_SIP4_and_PFIso40_and_PostTS", "DoubleMu17TkMu8_TkMu8_from_PF_and_SIP4_and_PFIso40_and_PostTS"
                                 ,"DoubleMu17Mu8_Mu17_EMU_from_PF_and_SIP4_and_PFIso40", "DoubleMu17Mu8_Mu17_EMU_from_PF_and_SIP4_and_PFIso40_and_PostTS"
                                 ,"HLT_OrMu17_from_PF_and_SIP4_and_PFIso40", "HLT_OrMu8_from_PF_and_SIP4_and_PFIso40"
                                 ,"HLT_Mu17_from_PF_and_SIP4_and_PFIso40_and_2011A", "HLT_Mu8_from_PF_and_SIP4_and_PFIso40_and_2011A"
                                 ,"HLT_Mu17_from_PF_and_SIP4_and_PFIso40_and_2011B", "HLT_Mu8_from_PF_and_SIP4_and_PFIso40_and_2011B"
                                 };
    const char *titles[nids] = { "PF ID", "S_{IP3D} < 4", "PF Rel Iso < 0.4"
                                 ,"PF ID", "PF ID"
                                 ,"Glb ID", "TM Arb ID", "ZZ Loose ID", "Tight ID (no IP)", "Tight ID (with IP)"
                                 ,"PF Rel Iso < 0.2 (#Delta#beta)", "PF Rel Iso < 0.12 (#Delta#beta)"
                                 ,"PF Rel Iso < 0.4", "PF Rel Iso < 0.2 (#Delta#beta)", ,"PF Rel Iso < 0.12 (#Delta#beta)"
                                 ,"HLT Mu17", "HLT Mu8" 
                                 ,"HLT Mu17", "HLT Mu8" 
                                 ,"HLT Mu17 Tk", "HLT TkMu8"
                                 ,"Post TS HLT Mu17", "Post TS HLT Mu8" 
                                 ,"Post TS HLT Mu17 Tk", "Post TS HLT TkMu8"
                                 ,"HLTEmu Mu17", "Post TS HLTEmu Mu17"
                                 ,"HLT Mu17 OR", "HLT Mu8 OR" 
                                 ,"HLT Mu17", "HLT Mu8" 
                                 ,"HLT Mu17", "HLT Mu8" 
                                 };
    for (int i = 0; i < nids; ++i) {
        TString idname(ids[i]);
        TString mcname = idname; 
        mcname.ReplaceAll("_and_PostTS",""); 
        mcname.ReplaceAll("_from_2011A",""); 
        mcname.ReplaceAll("_from_2011B",""); 
        mcname.ReplaceAll("_and_2011A",""); 
        mcname.ReplaceAll("_and_2011B",""); 
        retitle = TString(titles[i])+" efficiency";

        TDirectory *fit_pt_eta = gFile->GetDirectory(basedir+"/"+idname+"_pt_abseta/");
        if (fit_pt_eta) {
            yMin = 0.0; yMax = 1.1;
            if (ref != 0) {
                TDirectory *ref_pt_eta = ref->GetDirectory(basedir+"/"+mcname+"_pt_abseta/");
                extraSpam = "        |#eta| < 1.2"; refstack(fit_pt_eta, ref_pt_eta, idname+"_pt_barrel",  "pt_PLOT_abseta_bin0");
                extraSpam = "  1.2 < |#eta| < 2.4"; refstack(fit_pt_eta, ref_pt_eta, idname+"_pt_endcaps", "pt_PLOT_abseta_bin1");
            } else {
                extraSpam = "        |#eta| < 1.2"; single(fit_pt_eta, idname+"_pt_barrel",  "pt_PLOT_abseta_bin0");
                extraSpam = "  1.2 < |#eta| < 2.4"; single(fit_pt_eta, idname+"_pt_endcaps", "pt_PLOT_abseta_bin1");
            }
        }
        fit_pt_eta = gFile->GetDirectory(basedir+"/"+idname+"_pt_abseta2/");
        if (fit_pt_eta) {
            yMin = 0.0; yMax = 1.1;
            if (ref != 0) {
                TDirectory *ref_pt_eta = ref->GetDirectory(basedir+"/"+mcname+"_pt_abseta2/");
                extraSpam = "        |#eta| < 1.2"; refstack(fit_pt_eta, ref_pt_eta, idname+"_pt_barrel2",  "pt_PLOT_abseta_bin0");
                extraSpam = "  1.2 < |#eta| < 2.4"; refstack(fit_pt_eta, ref_pt_eta, idname+"_pt_endcaps2", "pt_PLOT_abseta_bin1");
            } else {
                extraSpam = "        |#eta| < 0.4"; single(fit_pt_eta, idname+"_pt_barrel2",  "pt_PLOT_abseta_bin0");
                extraSpam = "  1.2 < |#eta| < 2.4"; single(fit_pt_eta, idname+"_pt_endcaps2", "pt_PLOT_abseta_bin1");
            }
        }

        TDirectory *fit_eta = gFile->GetDirectory(basedir+"/"+idname+"_eta/");
        if (fit_eta) {
            yMin = 0.89; yMax = 1.019;
            if (retitle.Contains("HLT")) { yMin = 0.0; yMax = 1.1; }
            if (idname.BeginsWith("Glb")) { yMin = 0.7; yMax = 1.049; }
            if (idname.BeginsWith("PFIso12DB")) { yMin = 0.7; yMax = 1.049; }
            if (idname.BeginsWith("Tight")) { yMin = 0.7; yMax = 1.049; }
            if (ref != 0) {
                extraSpam = "    p_{T} > 20 GeV"; retitleX = "muon #eta";
                TDirectory *ref_eta = ref->GetDirectory(basedir+"/"+mcname+"_eta/");
                refstack(fit_eta, ref_eta, idname+"_eta",  "eta_PLOT");
            } else {
                single(fit_eta, idname+"_eta",  "eta_PLOT");
            }
        }
        fit_eta = gFile->GetDirectory(basedir+"/"+idname+"_eta2/");
        if (fit_eta) {
            yMin = 0.89; yMax = 1.019;
            if (retitle.Contains("HLT")) { yMin = 0.0; yMax = 1.1; }
            if (idname.BeginsWith("Glb")) { yMin = 0.7; yMax = 1.049; }
            if (idname.BeginsWith("PFIso12DB")) { yMin = 0.7; yMax = 1.049; }
            if (idname.BeginsWith("Tight")) { yMin = 0.7; yMax = 1.049; }
            if (ref != 0) {
                extraSpam = "    p_{T} > 20 GeV"; retitleX = "muon #eta";
                TDirectory *ref_eta = ref->GetDirectory(basedir+"/"+mcname+"_eta2/");
                refstack(fit_eta, ref_eta, idname+"_eta2",  "eta_PLOT");
            } else {
                single(fit_eta, idname+"_eta2",  "eta_PLOT");
            }
        }


        TDirectory *fit_vtx = gFile->GetDirectory(basedir+"/"+idname+"_vtx/");
        if (fit_vtx) {
            yMin = 0.89; yMax = 1.019;
            if (retitle.Contains("HLT")) { yMin = 0.0; yMax = 1.1; }
            if (idname.BeginsWith("PFIso12DB")) { yMin = 0.7; yMax = 1.049; }
            if (idname.BeginsWith("Tight")) { yMin = 0.7; yMax = 1.049; }
            TDirectory *ref_vtx = ref ? ref->GetDirectory(basedir+"/"+mcname+"_vtx/") : 0;
            retitleX = "Number of vertices";
            if (ref_vtx) {
                refstack(fit_vtx, ref_vtx, idname+"_vtx",  "tag_nVertices_PLOT");
            } else {
                single(fit_vtx, idname+"_vtx",  "tag_nVertices_PLOT");
            }
        }

        TDirectory *fit_overall_abseta = gFile->GetDirectory(basedir+"/"+idname+"_overall_abseta/");
        if (fit_overall_abseta) {
            yMin = 0.85; yMax = 1.019;
            if (retitle.Contains("HLT")) { yMin = 0.0; yMax = 1.1; }
            TDirectory *ref_overall_abseta = ref ? ref->GetDirectory(basedir+"/"+idname+"_overall_abseta/") : 0;
            if (ref_overall_abseta) {
                extraSpam = "    p_{T} > 20 GeV"; retitleX = "muon |#eta|";
                refstack(fit_overall_abseta, ref_overall_abseta, idname+"_overall_abseta",  "abseta_PLOT_");
            } else {
                single(fit_overall_abseta, idname+"_overall_abseta",  "abseta_PLOT_");
            }
        }
    }
}


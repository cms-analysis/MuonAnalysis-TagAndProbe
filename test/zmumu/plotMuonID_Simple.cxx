void plotMuonID_Simple() {
    // Open files
    TFile *fitIso = TFile::Open("TnP_Muon_Iso_Simple.root");
    TFile *fitID  = TFile::Open("TnP_Muon_ID_Simple.root");

    // Plot one of the canvases
    fitIso->Get("tpTree/Iso_vtx_tight/fit_eff_plots/tag_nVertices_PLOT_PF_pass_&_tag_IsoMu24_eta2p1_pass")->Draw();

    // Get a 1D efficiency plot as a TGraph
    TGraphAsymmErrors *eff_barrel  = (TGraphAsymmErrors *) fitID->Get("tpTree/PF_pt_eta/fit_eff_plots/pt_PLOT_abseta_bin0_&_tag_IsoMu24_eta2p1_pass")->FindObject("hxy_fit_eff");
    TGraphAsymmErrors *eff_endcaps = (TGraphAsymmErrors *) fitID->Get("tpTree/PF_pt_eta/fit_eff_plots/pt_PLOT_abseta_bin1_&_tag_IsoMu24_eta2p1_pass")->FindObject("hxy_fit_eff");
    TCanvas *c1 = new TCanvas("c1","c1"); c1->Divide(2,1);
    c1->cd(1); eff_barrel->Draw("AP");
    c1->cd(2); eff_endcaps->Draw("AP");

   // Print an efficiency table
    RooDataSet *dataset = (RooDataSet *) fitID->Get("tpTree/PF_pt_eta/fit_eff");
    for (int i = 0; i < dataset->numEntries(); ++i) {
        const RooArgSet &point = *dataset->get(i);
        RooRealVar &eta = point["abseta"], &pt = point["pt"], &eff = point["efficiency"];
        printf("bin %3d: abseta [%+5.3f, %+5.3f], pt [%5.1f, %5.1f]: efficiency: %6.2f %+5.2f/%+5.2f %%\n",
            i,
            eta.getVal() + eta.getAsymErrorLo(), eta.getVal() + eta.getAsymErrorHi(),
            pt.getVal() + pt.getAsymErrorLo(), pt.getVal() + pt.getAsymErrorHi(),
            100*eff.getVal(), 100*eff.getAsymErrorLo(), 100*eff.getAsymErrorHi());
    }

}

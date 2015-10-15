#include "TTree.h"
#include "TFile.h"
#include "TStopwatch.h"
#include "MuonEffectiveArea.h"

void addEAMiniIso() {
    TTree *tIn  = (TTree *) gFile->Get("tpTree/fitter_tree");
    Float_t pt, eta, chHad, nHad, phot, rho;
    tIn->SetBranchAddress("pt", &pt);
    tIn->SetBranchAddress("eta", &eta);
    tIn->SetBranchAddress("miniIsoCharged", &chHad);
    tIn->SetBranchAddress("miniIsoNeutrals", &nHad);
    tIn->SetBranchAddress("miniIsoPhotons",     &phot);
    tIn->SetBranchAddress("fixedGridRhoFastjetCentralNeutral",     &rho);

    TFile *fOut = new TFile("tnpZ_withEAMiniIso.root", "RECREATE");
    fOut->mkdir("tpTree")->cd();
    TTree *tOut = tIn->CloneTree(0);
    Float_t pfCombRelMiniIsoEACorr,pfCombAbsMiniIsoEACorr;
    tOut->Branch("pfCombAbsMiniIsoEACorr", &pfCombAbsMiniIsoEACorr, "pfCombAbsMiniIsoEACorr/F");
    tOut->Branch("pfCombRelMiniIsoEACorr", &pfCombRelMiniIsoEACorr, "pfCombRelMiniIsoEACorr/F");

    MuonEffectiveArea::MuonEffectiveAreaTarget effAreaTarget = MuonEffectiveArea::kMuEASpring15_25ns; // new 2015
    MuonEffectiveArea::MuonEffectiveAreaType   effAreaType   = MuonEffectiveArea::kMuMiniIso03;

    int step = tIn->GetEntries()/1000;
    double evDenom = 100.0/double(tIn->GetEntries());
    TStopwatch timer; timer.Start();
    for (int i = 0, n = tIn->GetEntries(); i < n; ++i) {
        tIn->GetEntry(i);
        Float_t ea_tot = MuonEffectiveArea::GetMuonEffectiveArea(effAreaType, fabs(eta), effAreaTarget);
        
        pfCombAbsMiniIsoEACorr = (chHad + max(0.0, nHad - rho * ea_tot * ((10.0/min(max(pt, 50),200))/0.3)**2)
        pfCombRelMiniIsoEACorr = pfCombAbsMiniIsoEACorr/pt;
        if (i < 20) {
            printf("muon with pt = %.2f, eta = %+5.2f:", pt, eta);
            printf("   charged hadrons %6.3f, neutral hadrons %6.3f, photons %6.3f ", chHad, nHad, phot);
            printf("   rho %6.3f, ea %6.3f", rho, ea_tot);
            printf("   pfCombAbsIsoEAcorr %6.3f\n", pfCombAbsMiniIsoEACorr);
            printf("   pfCombRelIsoEAcorr %6.3f\n", pfCombRelMiniIsoEACorr);
        }
        tOut->Fill();
        //if (i > 10000) break;
        if ((i+1) % step == 0) { 
            double totalTime = timer.RealTime()/60.; timer.Continue();
            double fraction = double(i+1)/double(n+1), remaining = totalTime*(1-fraction)/fraction;
            printf("Done %9d/%9d   %5.1f%%   (elapsed %5.1f min, remaining %5.1f min)\n", i, n, i*evDenom, totalTime, remaining); 
            fflush(stdout); 
        }
    }

    tOut->AutoSave(); // according to root tutorial this is the right thing to do
    fOut->Close();
}

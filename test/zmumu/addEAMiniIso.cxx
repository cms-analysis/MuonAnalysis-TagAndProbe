#include "TTree.h"
#include "TFile.h"
#include "TStopwatch.h"
#include "MuonEffectiveArea.h"

void addEAMiniIso() {
    TTree *tIn  = (TTree *) gFile->Get("tpTree/fitter_tree");
    Float_t pt, eta, chHad, nHad, phot, rho, activity_chHad, activity_PUchHad, activity_nHad, activity_phot;
    tIn->SetBranchAddress("pt", &pt);
    tIn->SetBranchAddress("eta", &eta);
    tIn->SetBranchAddress("miniIsoCharged", &chHad);
    tIn->SetBranchAddress("miniIsoNeutrals", &nHad);
    tIn->SetBranchAddress("miniIsoPhotons",     &phot);
    tIn->SetBranchAddress("activity_miniIsoCharged", &activity_chHad);
    tIn->SetBranchAddress("activity_miniIsoPUCharged", &activity_PUchHad);
    tIn->SetBranchAddress("activity_miniIsoNeutrals", &activity_nHad);
    tIn->SetBranchAddress("activity_miniIsoPhotons",     &activity_phot);
    tIn->SetBranchAddress("fixedGridRhoFastjetCentralNeutral",     &rho);

    TFile *fOut = new TFile("tnpZ_withEAMiniIso.root", "RECREATE");
    fOut->mkdir("tpTree")->cd();
    TTree *tOut = tIn->CloneTree(0);
    Float_t pfCombRelMiniIsoEACorr,pfCombAbsMiniIsoEACorr,pfCombRelActivitydBCorr,pfCombAbsActivitydBCorr;
    tOut->Branch("pfCombAbsMiniIsoEACorr", &pfCombAbsMiniIsoEACorr, "pfCombAbsMiniIsoEACorr/F");
    tOut->Branch("pfCombRelMiniIsoEACorr", &pfCombRelMiniIsoEACorr, "pfCombRelMiniIsoEACorr/F");
    tOut->Branch("pfCombRelActivitydBCorr", &pfCombRelActivitydBCorr, "pfCombRelActivitydBCorr/F");
    tOut->Branch("pfCombAbsActivitydBCorr", &pfCombAbsActivitydBCorr, "pfCombAbsActivitydBCorr/F");

    MuonEffectiveArea::MuonEffectiveAreaTarget effAreaTarget = MuonEffectiveArea::kMuEASpring15_25ns; // new 2015
    MuonEffectiveArea::MuonEffectiveAreaType   effAreaType   = MuonEffectiveArea::kMuMiniIso03;

    int step = tIn->GetEntries()/1000;
    double evDenom = 100.0/double(tIn->GetEntries());
    TStopwatch timer; timer.Start();
    for (int i = 0, n = tIn->GetEntries(); i < n; ++i) {
        tIn->GetEntry(i);
        Float_t ea_tot = MuonEffectiveArea::GetMuonEffectiveArea(effAreaType, fabs(eta), effAreaTarget);
        
        pfCombAbsMiniIsoEACorr = (chHad + max(0.0, nHad - rho * ea_tot * pow((10.0/min(max((double) pt, 50.),200.))/0.3,2)));
        pfCombRelMiniIsoEACorr = pfCombAbsMiniIsoEACorr/pt;
        
        pfCombAbsActivitydBCorr = (activity_chHad + max( (double) (activity_nHad +  activity_phot -  activity_PUchHad/2) , 0.0));
        pfCombRelActivitydBCorr = pfCombAbsActivitydBCorr/pt;
                
        if (i < 20) {
            printf("muon with pt = %.2f, eta = %+5.2f:", pt, eta);
            printf("   charged hadrons %6.3f, neutral hadrons %6.3f, photons %6.3f ", chHad, nHad, phot);
            printf("   rho %6.3f, ea %6.3f", rho, ea_tot);
            printf("   pfCombAbsMiniIsoEAcorr %6.3f\n", pfCombAbsMiniIsoEACorr);
            printf("   pfCombRelMiniIsoEAcorr %6.3f\n", pfCombRelMiniIsoEACorr);
            printf("   activity charged hadrons %6.3f, PU charged hadrons %6.3f, neutral hadrons %6.3f, photons %6.3f ", activity_chHad, activity_PUchHad, activity_nHad, activity_phot);
            printf("   pfCombAbsActivitydBCorr %6.3f\n", pfCombAbsActivitydBCorr);
            printf("   pfCombRelActivitydBCorr %6.3f\n", pfCombRelActivitydBCorr);
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

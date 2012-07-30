#include "TTree.h"
#include "TFile.h"
#include "TStopwatch.h"
#include "Muon/MuonAnalysisTools/interface/MuonEffectiveArea.h"

void addEAIso() {
    TTree *tIn  = (TTree *) gFile->Get("tpTree/fitter_tree");
    Float_t pt, eta, chHad, nHad, phot, rho;
    tIn->SetBranchAddress("pt", &pt);
    tIn->SetBranchAddress("eta", &eta);
    tIn->SetBranchAddress("chargedHadIso04", &chHad);
    tIn->SetBranchAddress("neutralHadIso04", &nHad);
    tIn->SetBranchAddress("photonIso04",     &phot);
    tIn->SetBranchAddress("kt6RhoNeu05",     &rho); // 2012
    //tIn->SetBranchAddress("kt6RhoAll",     &rho); // 2011

    TFile *fOut = new TFile("tnpZ_withEAIso.root", "RECREATE");
    fOut->mkdir("tpTree")->cd();
    TTree *tOut = tIn->CloneTree(0);
    Float_t pfCombRelIso04EACorr;
    tOut->Branch("pfCombRelIso04EACorr", &pfCombRelIso04EACorr, "pfCombRelIso04EACorr/F");

    //MuonEffectiveArea::MuonEffectiveAreaTarget effAreaTarget = MuonEffectiveArea::kMuEAData2011; 
    MuonEffectiveArea::MuonEffectiveAreaTarget effAreaTarget = MuonEffectiveArea::kMuEAData2012;
    MuonEffectiveArea::MuonEffectiveAreaType   effAreaType   = MuonEffectiveArea::kMuGammaAndNeutralHadronIso04;

    int step = tIn->GetEntries()/1000;
    double evDenom = 100.0/double(tIn->GetEntries());
    TStopwatch timer; timer.Start();
    for (int i = 0, n = tIn->GetEntries(); i < n; ++i) {
        tIn->GetEntry(i);
        Float_t ea_tot = MuonEffectiveArea::GetMuonEffectiveArea(effAreaType, fabs(eta), effAreaTarget);
        pfCombRelIso04EACorr = (chHad + max(0.f, nHad + phot - ea_tot*rho))/pt;
        if (i < 20) {
            printf("muon with pt = %.2f, eta = %+5.2f:", pt, eta);
            printf("   charged hadrons %6.3f, neutral hadrons %6.3f, photons %6.3f ", chHad, nHad, phot);
            printf("   rho %6.3f, ea %6.3f", rho, ea_tot);
            printf("   pfCombRelIsoEAcorr %6.3f\n", pfCombRelIso04EACorr);
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

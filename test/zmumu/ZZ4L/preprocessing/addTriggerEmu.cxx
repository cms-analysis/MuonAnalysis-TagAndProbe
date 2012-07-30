#include "TTree.h"
#include "TFile.h"
#include "TStopwatch.h"
#include "Muon/MuonAnalysisTools/interface/MuonEffectiveArea.h"

void addTriggerEmu() {
    TTree *tIn  = (TTree *) gFile->Get("tpTree/fitter_tree");
    Float_t l3pt; 
    Int_t DoubleMu17Mu8_Mu17, DoubleMu17Mu8_Mu8;
    Int_t DoubleMu17TkMu8_Mu17, DoubleMu17TkMu8_TkMu8;
    tIn->SetBranchAddress("l3pt", &l3pt);
    tIn->SetBranchAddress("DoubleMu17Mu8_Mu17", &DoubleMu17Mu8_Mu17);
    tIn->SetBranchAddress("DoubleMu17Mu8_Mu8", &DoubleMu17Mu8_Mu8);
    tIn->SetBranchAddress("DoubleMu17TkMu8_Mu17", &DoubleMu17TkMu8_Mu17);
    tIn->SetBranchAddress("DoubleMu17TkMu8_TkMu8", &DoubleMu17TkMu8_TkMu8);

    TFile *fOut = new TFile("tnpZ_withTriggerEmu.root", "RECREATE");
    fOut->mkdir("tpTree")->cd();
    TTree *tOut = tIn->CloneTree(0);
    Int_t Mu8, Mu17, TkMu8, TkMu17, OrMu8, OrMu17;
    tOut->Branch("HLT_Mu8", &Mu8, "HLT_Mu8/I");
    tOut->Branch("HLT_Mu17", &Mu17, "HLT_Mu17/I");
    tOut->Branch("HLT_TkMu8", &TkMu8, "HLT_TkMu8/I");
    tOut->Branch("HLT_TkMu17", &TkMu17, "HLT_TkMu17/I");
    tOut->Branch("HLT_OrMu8", &OrMu8, "HLT_OrMu8/I");
    tOut->Branch("HLT_OrMu17", &OrMu17, "HLT_OrMu17/I");

    int step = tIn->GetEntries()/1000;
    double evDenom = 100.0/double(tIn->GetEntries());
    TStopwatch timer; timer.Start();
    for (int i = 0, n = tIn->GetEntries(); i < n; ++i) {
        tIn->GetEntry(i);
        Mu17 = DoubleMu17Mu8_Mu17;
        Mu8  = DoubleMu17Mu8_Mu8;
        TkMu17 = DoubleMu17TkMu8_Mu17;
        TkMu8  = DoubleMu17TkMu8_TkMu8;
        OrMu17 = (Mu17 || TkMu17);
        OrMu8 = (Mu8 || TkMu8);
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

#include "TTree.h"
#include "TFile.h"
#include "TSystem.h"
#include "TStopwatch.h"
#include <string>
#include <vector>
#include <cstdlib>
#include <cmath>
#define STANDALONE
#include "Muon/MuonAnalysisTools/interface/MuonEffectiveArea.h"
#include "Muon/MuonAnalysisTools/interface/MuonMVAEstimator.h"

void addMVAIso() {
    const char *base=getenv("CMSSW_BASE");
    std::string baseFolder(base);
    baseFolder += "/src/Muon/MuonAnalysisTools/data/";
    std::vector<string> manualCatNonTrigWeigths;
    manualCatNonTrigWeigths.push_back(baseFolder+"/MuonIsoMVA_sixie-BarrelPt5To10_V0_BDTG.weights.xml");
    manualCatNonTrigWeigths.push_back(baseFolder+"/MuonIsoMVA_sixie-EndcapPt5To10_V0_BDTG.weights.xml");
    manualCatNonTrigWeigths.push_back(baseFolder+"/MuonIsoMVA_sixie-BarrelPt10ToInf_V0_BDTG.weights.xml");
    manualCatNonTrigWeigths.push_back(baseFolder+"/MuonIsoMVA_sixie-EndcapPt10ToInf_V0_BDTG.weights.xml");
    manualCatNonTrigWeigths.push_back(baseFolder+"/MuonIsoMVA_sixie-Tracker_V0_BDTG.weights.xml");
    manualCatNonTrigWeigths.push_back(baseFolder+"/MuonIsoMVA_sixie-Global_V0_BDTG.weights.xml");

    MuonMVAEstimator *muMVANonTrig  = new MuonMVAEstimator();
    muMVANonTrig->initialize("MuonIso_BDTG_IsoRings",MuonMVAEstimator::kIsoRings,true,manualCatNonTrigWeigths);
    muMVANonTrig->SetPrintMVADebug(kFALSE);

    TTree *tIn  = (TTree *) ((TFile*)gROOT->GetListOfFiles()->At(0))->Get("tpTree/fitter_tree");
    Float_t pt, eta, rho; Int_t global, tracker;
    Float_t ChargedIso_DR0p0To0p1, NeutralHadronIso_DR0p0To0p1, GammaIso_DR0p0To0p1;
    Float_t ChargedIso_DR0p1To0p2, NeutralHadronIso_DR0p1To0p2, GammaIso_DR0p1To0p2;
    Float_t ChargedIso_DR0p2To0p3, NeutralHadronIso_DR0p2To0p3, GammaIso_DR0p2To0p3;
    Float_t ChargedIso_DR0p3To0p4, NeutralHadronIso_DR0p3To0p4, GammaIso_DR0p3To0p4;
    Float_t ChargedIso_DR0p4To0p5, NeutralHadronIso_DR0p4To0p5, GammaIso_DR0p4To0p5;
    tIn->SetBranchAddress("pt",  &pt);
    tIn->SetBranchAddress("eta", &eta);
    tIn->SetBranchAddress("Glb", &global);
    tIn->SetBranchAddress("TM",  &tracker);
    tIn->SetBranchAddress("kt6RhoAll",     &rho);
    tIn->SetBranchAddress("ChargedIso_DR0p0To0p1", &ChargedIso_DR0p0To0p1);
    tIn->SetBranchAddress("ChargedIso_DR0p1To0p2", &ChargedIso_DR0p1To0p2);
    tIn->SetBranchAddress("ChargedIso_DR0p2To0p3", &ChargedIso_DR0p2To0p3);
    tIn->SetBranchAddress("ChargedIso_DR0p3To0p4", &ChargedIso_DR0p3To0p4);
    tIn->SetBranchAddress("ChargedIso_DR0p4To0p5", &ChargedIso_DR0p4To0p5);
    tIn->SetBranchAddress("NeutralHadronIso_DR0p0To0p1", &NeutralHadronIso_DR0p0To0p1);
    tIn->SetBranchAddress("NeutralHadronIso_DR0p1To0p2", &NeutralHadronIso_DR0p1To0p2);
    tIn->SetBranchAddress("NeutralHadronIso_DR0p2To0p3", &NeutralHadronIso_DR0p2To0p3);
    tIn->SetBranchAddress("NeutralHadronIso_DR0p3To0p4", &NeutralHadronIso_DR0p3To0p4);
    tIn->SetBranchAddress("NeutralHadronIso_DR0p4To0p5", &NeutralHadronIso_DR0p4To0p5);
    tIn->SetBranchAddress("GammaIso_DR0p0To0p1", &GammaIso_DR0p0To0p1);
    tIn->SetBranchAddress("GammaIso_DR0p1To0p2", &GammaIso_DR0p1To0p2);
    tIn->SetBranchAddress("GammaIso_DR0p2To0p3", &GammaIso_DR0p2To0p3);
    tIn->SetBranchAddress("GammaIso_DR0p3To0p4", &GammaIso_DR0p3To0p4);
    tIn->SetBranchAddress("GammaIso_DR0p4To0p5", &GammaIso_DR0p4To0p5);

    TFile *fOut = new TFile("tnpZ_withMVAIso.root", "RECREATE");
    fOut->mkdir("tpTree")->cd();
    TTree *tOut = tIn->CloneTree(0);
    Float_t mvaIso; Int_t mvaIsoCut;
    tOut->Branch("mvaIso", &mvaIso, "mvaIso/F");
    tOut->Branch("mvaIsoCut", &mvaIsoCut, "mvaIsoCut/I");

    MuonEffectiveArea::MuonEffectiveAreaTarget effAreaTarget = MuonEffectiveArea::kMuEAData2011; // NOTE: for MVA we use 2011 also for 2012 for now.

    int step = tIn->GetEntries()/1000;
    double evDenom = 100.0/double(tIn->GetEntries());
    TStopwatch timer; timer.Start();
    for (int i = 0, n = tIn->GetEntries(); i < n; ++i) {
        tIn->GetEntry(i);
        mvaIso = muMVANonTrig->mvaValue_Iso(pt, eta, global, tracker, rho, effAreaTarget,
                                                ChargedIso_DR0p0To0p1,
                                                ChargedIso_DR0p1To0p2,
                                                ChargedIso_DR0p2To0p3,
                                                ChargedIso_DR0p3To0p4,
                                                ChargedIso_DR0p4To0p5,
                                                GammaIso_DR0p0To0p1,
                                                GammaIso_DR0p1To0p2,
                                                GammaIso_DR0p2To0p3,
                                                GammaIso_DR0p3To0p4,
                                                GammaIso_DR0p4To0p5,
                                                NeutralHadronIso_DR0p0To0p1,
                                                NeutralHadronIso_DR0p1To0p2,
                                                NeutralHadronIso_DR0p2To0p3,
                                                NeutralHadronIso_DR0p3To0p4,
                                                NeutralHadronIso_DR0p4To0p5,
                                                false);
        if (global && tracker) {
            if (pt <  10 && abs(eta) <  1.5)   mvaIsoCut = (mvaIso > -0.593);
            if (pt >= 10 && abs(eta) <  1.5)   mvaIsoCut = (mvaIso >  0.337);
            if (pt <  10 && abs(eta) >= 1.5)   mvaIsoCut = (mvaIso > -0.767);
            if (pt >= 10 && abs(eta) >= 1.5)   mvaIsoCut = (mvaIso >  0.410);
        } else if (!global && tracker) {
            mvaIsoCut = (mvaIso > -0.989);
        } else if (tracker && !global) {
            mvaIsoCut = (mvaIso > -0.995);
        } else {
            mvaIsoCut = 0;
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

    tOut->Write();
    fOut->Close();
}

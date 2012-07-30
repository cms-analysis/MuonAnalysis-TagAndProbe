void subTree(TString dir="tpTree", TString cut="(tag_IsoMu24_eta2p1 || tag_IsoMu24) && tag_chargedHadIso04/tag_pt < 0.2", TString newFile="subTree_IsoMu24.root") {
    TTree *in  = (TTree *)gFile->Get(dir+"/fitter_tree");
    TFile *fout = new TFile(newFile, "RECREATE");
    TDirectory *dout = fout->mkdir(dir); dout->cd();

    TTree *out = in->CopyTree(cut);
    std::cout << "INPUT TREE (" << in->GetEntries() << " ENTRIES)" << std::endl;
    //in->Print();
    std::cout << "OUTPUT TREE (" << out->GetEntries() << " ENTRIES)" << std::endl;
    //out->Print();
    dout->WriteTObject(out, "fitter_tree");
    fout->Close();
}

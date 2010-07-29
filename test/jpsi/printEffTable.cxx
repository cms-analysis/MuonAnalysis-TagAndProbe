void printEffTable(TString dir, TString meas, TString vars="pt:eta", TString method="fit") {
    RooDataSet *ds = (RooDataSet *) gFile->Get(dir+"/"+meas+"/"+method+"_eff");
    if (ds == 0) {
        std::cerr << "NOT FOUND: " << (dir+"/"+meas+"/"+method+"_eff") << std::endl;
        if (gFile->GetDirectory(dir) != 0) {
            gFile->GetDirectory(dir)->ls();
        } else {
            gFile->ls();
        }
    }
    ds->tree()->Scan("efficiency:efficiency_aerr_lo:efficiency_aerr_hi:"+vars);
}

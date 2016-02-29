#include "../plotUtil_Tracking.cxx"

void plotFits(TString FileName = "../fits/Run2015D_16Dec2015-v1_Json/TnP_Tracking_dr030e015_TestData_Run2015D_16Dec2015-v1_Json_four.root"){

  cmslabel = "CMS Preliminary";
  cmsconditionslabel = "#sqrt{s} = 13 TeV, L = 2.23 fb^{-1}";

  TFile miofile(FileName,"read");
  TCanvas *c1 = (TCanvas*)miofile.Get("tpTreeSta/eff_four_dr030e015/abseta_bin0__pt_bin0__MuX_L2Mu0_L2_pass__outerValidHits_pass__tag_Mu5_L2Mu0_Mu_pass__gaussPlusCubic/fit_canvas");
  TCanvas *c_pf = new TCanvas("c_pf");
  TObject *obj; 
  TIter next(c1->GetListOfPrimitives());
  while ((obj=next())) {
    cout << "Reading: "  << obj->ClassName() << " - "<< obj->GetName() << endl;
    if (obj->InheritsFrom("TPad")) {
      cout << "tpad: " << obj->GetName() << endl;
      TString name = obj->GetName();
//    if(name.Contains("_1")){
      cout << "Plotting: "  << obj->GetName() << endl;
      TCanvas *c_pp = new TCanvas(obj->GetName());
      c_pp->cd(); 
      //obj->Draw();
      TPad * p = obj->Clone();	//working on RooFit v3.60, not sure in other root version
      p->SetPad("pp", "pp",0.01,0.01,0.99,0.99);
      p->SetFillColor(0);
      p->SetBorderMode(0);
      p->SetBorderSize(2);
      p->SetFrameBorderMode(0);
      p->SetFrameBorderMode(0);
      //p->DeleteExec("title");
      TIter next2(p->GetListOfPrimitives());
      TObject *obj2; 
      while ((obj2=next2())) {
        cout << "Reading in tpad: " << obj2->ClassName() << " - " << obj2->GetName() << endl;
        if (obj2->InheritsFrom("TPaveText")) {
          std::cout << "deleting TPaveText\n";
          p->DeleteToolTip(obj2);
        }
//      else if (obj2->InheritsFrom("TH1D")) {
//        std::cout << "working in TH1\n";
//        TH1D * histo = (TH1D*)p->GetPrimitive("frame_7f8a3d043100");
//        histo->SetTitle("");
//      }
        p->Draw();
             
      }
//        printcmsprelim();
//        printcmsconditionslabel();
//    c_pp.Close();

    }
  } 
}

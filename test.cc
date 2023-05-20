#include"SKFlatPlotter.cc"
void test(){
  SKFlatPlotter aa("FakeAnalyzer","data-mi-tau_mi-vv-wjets-tttw-aa");

  TString region="cpt/noZ/";

  TH1* hnum=aa.GetHist(0,"mM2018/"+region+"ss_l0etapt","project:y widthweight");
  TH1* hden=aa.GetHist(0,"MM2018/"+region+"ss_l0etapt","project:y widthweight");
  hden->SetLineColor(2);
  new TCanvas;
  hnum->Clone()->Draw();
  hden->Draw("same");

  new TCanvas;
  hnum->Divide(hden);
  hnum->Draw("hist e");
  hnum->GetYaxis()->SetRangeUser(0,2);

  hnum=aa.GetHist(0,"Mm2018/"+region+"ss_l1etapt","project:y widthweight");
  hden=aa.GetHist(0,"MM2018/"+region+"ss_l1etapt","project:y widthweight");
  hnum->Divide(hden);
  hnum->Draw("same hist e");
  hnum->SetLineColor(2);
  hnum->SetMarkerColor(2);
  gPad->SetLogx();

}
  

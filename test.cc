#include"SKFlatPlotter.cc"
void test(){
  SKFlatPlotter aa("FakeAnalyzer","data-mi-tau_mi-vv-wjets-tttw-aa");

  TH1* hnum=aa.GetHist(0,"mM2018/ss_l0etae_noZ","project:y widthweight");
  TH1* hden=aa.GetHist(0,"MM2018/ss_al0etae_noZ","project:y widthweight");
  hden->SetLineColor(2);
  new TCanvas;
  hnum->Clone()->Draw();
  hden->Draw("same");

  new TCanvas;
  hnum->Divide(hden);
  hnum->Draw("hist e");
  hnum->GetYaxis()->SetRangeUser(0,2);

  hnum=aa.GetHist(0,"Mm2018/ss_l0etae_noZ","project:y widthweight");
  hden=aa.GetHist(0,"MM2018/ss_al1etae_noZ","project:y widthweight");
  hnum->Divide(hden);
  hnum->Draw("same hist e");
  hnum->SetLineColor(2);
  hnum->SetMarkerColor(2);
  gPad->SetLogx();

}
  

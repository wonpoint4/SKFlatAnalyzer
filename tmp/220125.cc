Aepcor aa;

double RocEle2018_data_pt(double x,double eta=0.0,double r9=1.0){
  return aa.kScaleDT(x,eta,0.1,r9,0);
}
void DrawRocEle(){
  aa.init("data/Run2UltraLegacy_v2/2018/RoccoR/e_18UL.txt");
  TF1* f=new TF1("f","RocEle2018_data_pt(x)",10,1000);
  f->Draw();
  TF1* f2=new TF1("f2","RocEle2018_data_pt(x,2.5)",10,1000);
  f2->SetLineColor(3);
  f2->Draw("same");
  TF1* f3=new TF1("f3","RocEle2018_data_pt(x,0,0.2)",10,1000);
  f3->SetLineColor(4);
  f3->Draw("same");
  gPad->SetLogx();
  TH1* hist=(TH1*)f->GetHistogram();
  hist->GetYaxis()->SetRangeUser(0.9,1.1);
  hist->GetYaxis()->SetTitle("Rochester electron correction factor");
  hist->GetXaxis()->SetTitle("Electron p_{T} [GeV]");
  TLegend *leg=new TLegend(0.89,0.89,0.6,0.6);
  leg->AddEntry(f,"#eta=0.0, R9=1.0");
  leg->AddEntry(f2,"#eta=2.4, R9=1.0");
  leg->AddEntry(f3,"#eta=0.0, R9=0.2");
  leg->Draw();
  
}

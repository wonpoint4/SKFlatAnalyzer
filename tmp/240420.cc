#include"Plotter.cc"
void Draw(double ap=0.75,double am=0.75){
  Plotter::SetupStyle();
  TCanvas* c=gROOT->MakeDefCanvas();
  c->SetFrameLineWidth(0);
  TH1D* hist=new TH1D("hist","",4,0,4);
  hist->SetStats(0);
  hist->GetXaxis()->SetBinLabel(1,"#minus #minus");
  hist->GetXaxis()->SetBinLabel(2,"#plus #minus");
  hist->GetXaxis()->SetBinLabel(3,"#minus #plus");
  hist->GetXaxis()->SetBinLabel(4,"#plus #plus");
  hist->GetXaxis()->SetTitle("bjet-pair charge");
  hist->GetXaxis()->SetLabelSize(hist->GetXaxis()->GetLabelSize()*2);
  hist->GetYaxis()->SetLabelSize(0);
  hist->GetYaxis()->SetRangeUser(0,0.6);
  hist->GetYaxis()->SetAxisColor(0);
  hist->SetLineColor(kBlue-7);
  hist->SetLineWidth(3);
  double mm=1*am*(1-ap);
  double pm=0.5*ap*am+0.5*(1-ap)*(1-am);
  double mp=0.5*ap*am+0.5*(1-ap)*(1-am);
  double pp=1*ap*(1-am);
  cout<<mm+mp+pm+pp<<endl;
  hist->SetBinContent(1,mm);
  hist->SetBinContent(2,pm);
  hist->SetBinContent(3,mp);
  hist->SetBinContent(4,pp);
  hist->Draw();
  return c;
}

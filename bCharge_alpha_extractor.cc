#include "canvas_margin.h"
#include "kinFitPlotter.cc"
#if __has_include("dybAnalyzer.h")
#include "dybAnalyzer.h"
#endif

TH1* hist;
kinFitPlotter a("data mc");
kinFitPlotter b("ttlj");
void Getalphas_0D(TString channel="[em]2017/", TString suffix="");
void Getalphas_1D(TString channel="[em]2017/", TString suffix="");
double Getfc(TString histName);
double GetNorm(TString histName);

void bCharge_alpha_extractor(bool is1D = true){
  TString channel = "[em]2017/";

  if(!is1D) Getalphas_0D(channel);
  else Getalphas_1D(channel);
}

double Getfc(TString histName){
  double correct = b.GetHist(0, histName)->Integral();
  double wrong = b.GetHist(1, histName)->Integral();

  return correct / (correct + wrong);
}
double GetNorm(TString histName){
  double data = a.GetHist(0, histName)->Integral();
  double mc = a.GetHist(1, histName)->Integral();

  return data / mc;
}

void Getalphas_0D(TString channel, TString suffix = ""){
  cout<<"Get alpha in 0D"<<endl;
  TString lepb_histName = channel+"lepbjetChargeEasy_Lm";
  TString hadb_histName = channel+"hadbjetChargeEasy_Lm";

  double fc = Getfc(lepb_histName);
  double norm = GetNorm(lepb_histName);
  cout<<"fc = "<<fc<<", norm = "<<norm<<"\n"<<endl;

  kinFitPlotter c("data_sub ttlj", norm);

  double fp_lep_data = c.GetHist(0, lepb_histName)->GetBinContent(2) / (c.GetHist(0, lepb_histName)->GetBinContent(1) + c.GetHist(0, lepb_histName)->GetBinContent(2));
  double fm_had_data = c.GetHist(0, hadb_histName)->GetBinContent(1) / (c.GetHist(0, hadb_histName)->GetBinContent(1) + c.GetHist(0, hadb_histName)->GetBinContent(2));
  double fp_lep_mc = c.GetHist(1, lepb_histName)->GetBinContent(2) / (c.GetHist(1, lepb_histName)->GetBinContent(1) + c.GetHist(1, lepb_histName)->GetBinContent(2));
  double fm_had_mc = c.GetHist(1, hadb_histName)->GetBinContent(1) / (c.GetHist(1, hadb_histName)->GetBinContent(1) + c.GetHist(1, hadb_histName)->GetBinContent(2));

  cout<<"fp_lep_data = "<<fp_lep_data<<", fm_had_data = "<<fm_had_data<<endl;
  cout<<"fp_lep_mc = "<<fp_lep_mc<<", fm_had_mc = "<<fm_had_mc<<"\n"<<endl;

  double ap_data = (fc * (fp_lep_data + fc -1) + (1 - fc) * (fm_had_data + fc -1)) / (2 * fc - 1);
  double am_data = ((1 - fc) * (fp_lep_data + fc -1) + fc * (fm_had_data + fc -1)) / (2 * fc - 1);
  double ap_mc = (fc * (fp_lep_mc + fc -1) + (1 - fc) * (fm_had_mc + fc -1)) / (2 * fc - 1);
  double am_mc = ((1 - fc) * (fp_lep_mc + fc -1) + fc * (fm_had_mc + fc -1)) / (2 * fc - 1);

  cout<<"ap_data = "<<ap_data<<", am_data = "<<am_data<<endl;
  cout<<"ap_mc = "<<ap_mc<<", am_mc = "<<am_mc<<"\n"<<endl;
}

void Getalphas_1D(TString channel, TString suffix = ""){
  cout<<"Get alpha in 1D"<<endl;

  double fc = Getfc(channel+"lepbjetChargeEasy_Lm");
  for(unsigned int ch=0; ch<dybAnalyzer::afb_chbinnum; ch++){
    cout<<"\n"<<ch<<"th Charge bin "<<endl;
    TString lepb_histName = Form(channel+"lepbjetCharge%dEasy_Lm", ch);
    TString hadb_histName = Form(channel+"hadbjetCharge%dEasy_Lm", ch);

    double fc_bin = Getfc(lepb_histName);
    double norm = GetNorm(lepb_histName);
    cout<<"fc used = "<<fc<<", fc_bin = "<<fc_bin<<", norm = "<<norm<<"\n"<<endl;

    kinFitPlotter c("data_sub ttlj", norm);

    double fp_lep_data = c.GetHist(0, lepb_histName)->GetBinContent(2) / (c.GetHist(0, lepb_histName)->GetBinContent(1) + c.GetHist(0, lepb_histName)->GetBinContent(2));
    double fm_had_data = c.GetHist(0, hadb_histName)->GetBinContent(1) / (c.GetHist(0, hadb_histName)->GetBinContent(1) + c.GetHist(0, hadb_histName)->GetBinContent(2));
    double fp_lep_mc = c.GetHist(1, lepb_histName)->GetBinContent(2) / (c.GetHist(1, lepb_histName)->GetBinContent(1) + c.GetHist(1, lepb_histName)->GetBinContent(2));
    double fm_had_mc = c.GetHist(1, hadb_histName)->GetBinContent(1) / (c.GetHist(1, hadb_histName)->GetBinContent(1) + c.GetHist(1, hadb_histName)->GetBinContent(2));

    cout<<"fp_lep_data = "<<fp_lep_data<<", fm_had_data = "<<fm_had_data<<endl;
    cout<<"fp_lep_mc = "<<fp_lep_mc<<", fm_had_mc = "<<fm_had_mc<<"\n"<<endl;

    double ap_data = (fc * (fp_lep_data + fc -1) + (1 - fc) * (fm_had_data + fc -1)) / (2 * fc - 1);
    double am_data = ((1 - fc) * (fp_lep_data + fc -1) + fc * (fm_had_data + fc -1)) / (2 * fc - 1);
    double ap_mc = (fc * (fp_lep_mc + fc -1) + (1 - fc) * (fm_had_mc + fc -1)) / (2 * fc - 1);
    double am_mc = ((1 - fc) * (fp_lep_mc + fc -1) + fc * (fm_had_mc + fc -1)) / (2 * fc - 1);

    cout<<"ap_data = "<<ap_data<<", am_data = "<<am_data<<endl;
    cout<<"ap_mc = "<<ap_mc<<", am_mc = "<<am_mc<<"\n"<<endl;
  }
}

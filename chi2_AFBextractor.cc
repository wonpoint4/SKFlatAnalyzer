#include "canvas_margin.h"
#include "dybPlotter.cc"

dybPlotter a;
TH1* AFB;
vector<double> sth2_values = {0.23151, 0.23154, 0.23157, 0.2230, 0.2300, 0.2305, 0.2310, 0.2315, 0.2320, 0.2325, 0.2330 };
void Hists_1D_AFB(TString channel="[em][em]201[6-8][ab]?/m[52,200]/y[0,5]/", TString suffix="");
void Hists_2D_AFB(TString channel="[em][em]201[6-8][ab]?/m[52,200]/", TString frame="AFBrecoil");
void Plots_1D_chi2(TString inputfile="1D_sintheta.root", TString suffix="");
void Plots_2D_chi2(TString inputfile="2D_sintheta.root");
void Charge_Purity_Calc(TString channel="[em][em]201[6-8][ab]?/");
void Charge_Purity_Calc2(TString channel="[em][em]201[6-8][ab]?/");

void chi2_AFBextractor(){
  vector<TString> suffixes = {"_Tight1b_alljets", "_Tight1b_lepvetojets", "_Medium1b", "_Tightnb", "_Tight1b", "_Veto2b", "_Veto2j", "_MET75", "_ZbdPhi1p6", "_ZbpT60", ""};
  //for(const auto& suffix:suffixes){
  //  Hists_1D_AFB(suffix, "[em][em]201[6-8][ab]?/m[52,150]/y[0.2,5]/");
  //  Plots_1D_chi2(suffix);
  //}

  //Hists_1D_AFB();
  //Hists_2D_AFB();

  Plots_1D_chi2();
  //Plots_2D_chi2();

  //Charge_Purity_Calc();
  //Charge_Purity_Calc2();
}

void Hists_1D_AFB(TString channel, TString suffix){

  AFB = a.GetHist(0, channel+"AFBrecoil"+suffix+"(x)", "");
  AFB->SetName("Data_sthw2"+suffix);
  AFB->SaveAs("sthw2_Data"+suffix+".root");
  AFB = a.GetHist(1, channel+"AFBrecoil"+suffix+"(x)", "");
  AFB->SetName("MC_sthw2"+suffix);
  AFB->SaveAs("sthw2_MC"+suffix+".root");

  for(unsigned int i=0; i<sth2_values.size(); i++){
    AFB = a.GetHist(1, channel+"AFBrecoil"+suffix+"(x)", Form("suffix:_sthw2_%d:dy",i));
    AFB->SetName(Form("MC_sthw2"+suffix+"_%d",i));
    AFB->SaveAs(Form("sthw2_MC"+suffix+"_%d.root",i));
  }
}

vector<TString> chargebins = {"0", "0.1", "0.2", "0.6", "1", "3", "5"}; // should be equal with dybAnalyzer::afb_chbin

void Hists_2D_AFB(TString channel="[em][em]201[6-8][ab]?/m[52,200]/", TString frame="AFBrecoil"){
  for(unsigned int ch=0; ch<chargebins.size()-1; ch++){
    //cout<<channel+"y["+chargebins.at(ch)+","+chargebins.at(ch+1)+"]/"+frame+"(x)"<<endl;
    AFB = a.GetHist(0, channel+"y["+chargebins.at(ch)+","+chargebins.at(ch+1)+"]/"+frame+"(x)", "");
    AFB->SetName(Form("Data_sthw2_ch%d",ch));
    AFB->SaveAs(Form("sthw2_ch%d_Data.root",ch));
    AFB = a.GetHist(1, channel+"y["+chargebins.at(ch)+","+chargebins.at(ch+1)+"]/"+frame+"(x)", "");
    AFB->SetName(Form("MC_sthw2_ch%d",ch));
    AFB->SaveAs(Form("sthw2_ch%d_MC.root",ch));

    for(unsigned int i=0; i<sth2_values.size(); i++){
      AFB = a.GetHist(1, channel+"y["+chargebins.at(ch)+","+chargebins.at(ch+1)+"]/"+frame+"(x)", Form("suffix:_sthw2_%d:dy",i));
      AFB->SetName(Form("MC_sthw2_%d_ch%d",i,ch));
      AFB->SaveAs(Form("sthw2_%d_ch%d_MC.root",i,ch));
    }
  }
}

void Plots_1D_chi2(TString inputfile, TString suffix){

  TFile *file_AFB = new TFile(inputfile);

  TCanvas *c_AFB = new TCanvas("c_AFB", "", 1000, 1000);
  c_AFB->Draw();

  TH1D *hist_data = (TH1D*)file_AFB->Get("Data_sthw2"+suffix);
  TH1D *hist_mc   = (TH1D*)file_AFB->Get("MC_sthw2"+suffix);

  hist_mc->Draw();
  hist_data->Draw("same");

  TLegend *lg = new TLegend(0.6, 0.3, 0.9, 0.5);
  lg->SetFillStyle(0);
  lg->SetBorderSize(0);
  lg->AddEntry(hist_data, "Data set to MC", "lp");
  lg->AddEntry(hist_mc, "MC Nominal", "lp");
  lg->Draw();

  vector<double> chi2_values = {};

  for(unsigned int i=0; i<sth2_values.size(); i++){
    TH1D *hist_mc_sthw2 = (TH1D*)file_AFB->Get(Form("MC_sthw2"+suffix+"_%d",i));
    double chi2_tot = 0;

    for(unsigned int j=1; j<hist_mc_sthw2->GetNbinsX()-1; j++){ // GetNbinsX()+1
      hist_data->SetBinContent(j,hist_mc->GetBinContent(j)); // Use MC AFB values due to blinding

      double diff = hist_mc_sthw2->GetBinContent(j) - hist_data->GetBinContent(j);
      double sigma = sqrt(hist_data->GetBinError(j)*hist_data->GetBinError(j) + hist_mc_sthw2->GetBinError(j)*hist_mc_sthw2->GetBinError(j));
      double chi2 = diff*diff/(sigma*sigma);
      //cout<<diff<<", "<<sigma<<", "<<chi2<<endl;

      chi2_tot += chi2;
      hist_mc_sthw2->SetBinError(j,1e-10);
      hist_mc_sthw2->SetLineWidth(2);
    }
    cout<<"sthw2 "<<i<<"th("<<sth2_values.at(i)<<") chi2 = "<<chi2_tot<<endl;
    chi2_values.push_back(chi2_tot);

    if(i==3){
      hist_mc_sthw2->SetLineColor(kGreen+2);
      hist_mc_sthw2->Draw("same");
      lg->AddEntry(hist_mc_sthw2, "MC sin^{2}#theta = 0.22300", "lp");
    }else if(i==10){
      hist_mc_sthw2->SetLineColor(kBlue);
      hist_mc_sthw2->Draw("same");
      lg->AddEntry(hist_mc_sthw2, "MC sin^{2}#theta = 0.23300", "lp");
    }
  }

  // To print chi2 in easy format (python list)
  vector<unsigned int> chi2_index = {3, 4, 5, 6, 7, 0, 1, 2, 8, 9, 10};
  cout<<"\nchi2 = np.array([";
  for(unsigned int i=0; i<chi2_index.size(); i++){
    cout<<chi2_values.at(chi2_index.at(i))<<", ";
  }
  cout<<"])"<<endl;
  c_AFB->SaveAs("1D_AFBrecoil_plots"+suffix+".png");
}

void Plots_2D_chi2(TString inputfile){

  TFile *file_AFB = new TFile(inputfile);

  for(unsigned int ch=0; ch<chargebins.size()-1; ch++){
    TCanvas *c_AFB = new TCanvas("c_AFB", "", 1000, 1000);
    c_AFB->Draw();

    TH1D *hist_data = (TH1D*)file_AFB->Get(Form("Data_sthw2_ch%d",ch));
    TH1D *hist_mc   = (TH1D*)file_AFB->Get(Form("MC_sthw2_ch%d",ch));

    hist_mc->Draw();
    hist_data->Draw("same");

    TLegend *lg = new TLegend(0.6, 0.3, 0.9, 0.5);
    lg->SetFillStyle(0);
    lg->SetBorderSize(0);
    lg->AddEntry(hist_data, "Data set to MC", "lp");
    lg->AddEntry(hist_mc, "MC Nominal", "lp");
    lg->Draw();

    vector<double> chi2_values = {};

    for(unsigned int i=0; i<sth2_values.size(); i++){
      TH1D *hist_mc_sthw2 = (TH1D*)file_AFB->Get(Form("MC_sthw2_%d_ch%d",i,ch));
      double chi2_tot = 0;

      for(unsigned int j=1; j<hist_mc_sthw2->GetNbinsX()-1; j++){ // GetNbinsX()+1
        hist_data->SetBinContent(j,hist_mc->GetBinContent(j)); // Use MC AFB values due to blinding

	double diff = hist_mc_sthw2->GetBinContent(j) - hist_data->GetBinContent(j);
        double sigma = sqrt(hist_data->GetBinError(j)*hist_data->GetBinError(j) + hist_mc_sthw2->GetBinError(j)*hist_mc_sthw2->GetBinError(j));
        double chi2 = diff*diff/(sigma*sigma);
        //cout<<diff<<", "<<sigma<<", "<<chi2<<endl;

        chi2_tot += chi2;
        hist_mc_sthw2->SetBinError(j,1e-10);
        hist_mc_sthw2->SetLineWidth(2);
      }
      cout<<"sthw2 "<<i<<"th("<<sth2_values.at(i)<<") chi2 = "<<chi2_tot<<endl;
      chi2_values.push_back(chi2_tot);

      if(i==3){
        hist_mc_sthw2->SetLineColor(kGreen+2);
        hist_mc_sthw2->Draw("same");
        lg->AddEntry(hist_mc_sthw2, "MC sin^{2}#theta = 0.22300", "lp");
      }else if(i==10){
        hist_mc_sthw2->SetLineColor(kBlue);
        hist_mc_sthw2->Draw("same");
        lg->AddEntry(hist_mc_sthw2, "MC sin^{2}#theta = 0.23300", "lp");
      }
    }

    // To print chi2 in easy format (python list)
    vector<unsigned int> chi2_index = {3, 4, 5, 6, 7, 0, 1, 2, 8, 9, 10};
    cout<<"\nchi2 = np.array([";
    for(unsigned int i=0; i<chi2_index.size(); i++){
      cout<<chi2_values.at(chi2_index.at(i))<<", ";
    }
    cout<<"])"<<endl;
    c_AFB->SaveAs(Form("2D_AFBrecoil_plots_ch%d.png",ch));
  }
}

void Charge_Purity_Calc(TString channel="[em][em]bx201[6-8][ab]?/"){
  dybPlotter b("Data ^Dyb_mi+Dybbar_mi+mi+Dyc_mi+Dycbar_mi+Dyudsg_mi+tau_mi+vv+ss_mi+aa+tttw Dyb_mi Dybbar_mi");
  vector<TString> charges = {"P", "M"};

  for(unsigned int ch=0; ch<charges.size(); ch++){
    TString hist_name = channel+"bjetCharge_"+charges.at(ch)+"[2-5]";
    TH1* Data_hist  = b.GetHist(0, hist_name, "");
    TH1* allMC_hist = b.GetTH1(b.GetHist(1, hist_name, ""));
    TH1* Dyb_hist   = b.GetHist(2, hist_name, "");
    TH1* Dybbar_hist= b.GetHist(3, hist_name, "");

    double Data  = Data_hist->Integral();
    double allMC = allMC_hist->Integral();
    double Dyb   = Dyb_hist->Integral();
    double Dybbar= Dybbar_hist->Integral();

    cout<<"In "+hist_name+", Data = "<<Data<<", All MC = "<<allMC<<", DY+b = "<<Dyb<<", DY+bbar = "<<Dybbar<<", and Charge = "<<(charges.at(ch)=="P"? Dybbar: Dyb)/(Dyb+Dybbar)<<", and Purity = "<<(Dyb+Dybbar)/allMC<<", and Norm = "<<Data/allMC<<endl;

    for(unsigned int i=0; i<6; i++){
      hist_name = channel+"bjetCharge_"+charges.at(ch)+Form("%d",i);
      Data_hist  = b.GetHist(0, hist_name, "");
      allMC_hist = b.GetTH1(b.GetHist(1, hist_name, ""));
      Dyb_hist   = b.GetHist(2, hist_name, "");
      Dybbar_hist= b.GetHist(3, hist_name, "");

      Data  = Data_hist->Integral();
      allMC = allMC_hist->Integral();
      Dyb   = Dyb_hist->Integral();
      Dybbar= Dybbar_hist->Integral();

      cout<<"In "+hist_name+", Data = "<<Data<<", All MC = "<<allMC<<", DY+b = "<<Dyb<<", DY+bbar = "<<Dybbar<<", and Charge = "<<(charges.at(ch)=="P"? Dybbar: Dyb)/(Dyb+Dybbar)<<", and Purity = "<<(Dyb+Dybbar)/allMC<<", and Norm = "<<Data/allMC<<endl;
    }
  }
}

void Charge_Purity_Calc2(TString channel="[em][em]bx201[6-8][ab]?/"){
  dybPlotter b("Data ^Dyb_mi+Dybbar_mi+mi+Dyc_mi+Dycbar_mi+Dyudsg_mi+tau_mi+vv+ss_mi+aa+tttw Dyb_mi Dybbar_mi");
  vector<TString> charges = {"P", "M"};
  vector<TString> pts = {"L", "M", "H"};
  vector<TString> etas = {"B", "E"};

  // Pt Bins
  for(unsigned int pt=0; pt<pts.size(); pt++){
    for(unsigned int ch=0; ch<charges.size(); ch++){
      TString hist_name = channel+"bjetCharge_"+charges.at(ch)+pts.at(pt)+"[BE]";
      TH1* Data_hist  = b.GetHist(0, hist_name, "");
      TH1* allMC_hist = b.GetTH1(b.GetHist(1, hist_name, ""));
      TH1* Dyb_hist   = b.GetHist(2, hist_name, "");
      TH1* Dybbar_hist= b.GetHist(3, hist_name, "");

      double Data  = Data_hist->Integral();
      double allMC = allMC_hist->Integral();
      double Dyb   = Dyb_hist->Integral();
      double Dybbar= Dybbar_hist->Integral();

      cout<<"In "+hist_name+", Data = "<<Data<<", All MC = "<<allMC<<", DY+b = "<<Dyb<<", DY+bbar = "<<Dybbar<<", and Charge = "<<(charges.at(ch)=="P"? Dybbar: Dyb)/(Dyb+Dybbar)<<", and Purity = "<<(Dyb+Dybbar)/allMC<<", and Norm = "<<Data/allMC<<endl;
    }
  }

  // Eta Bins
  for(unsigned int eta=0; eta<etas.size(); eta++){
    for(unsigned int ch=0; ch<charges.size(); ch++){
      TString hist_name = channel+"bjetCharge_"+charges.at(ch)+"[LMH]"+etas.at(eta);
      TH1* Data_hist  = b.GetHist(0, hist_name, "");
      TH1* allMC_hist = b.GetTH1(b.GetHist(1, hist_name, ""));
      TH1* Dyb_hist   = b.GetHist(2, hist_name, "");
      TH1* Dybbar_hist= b.GetHist(3, hist_name, "");

      double Data  = Data_hist->Integral();
      double allMC = allMC_hist->Integral();
      double Dyb   = Dyb_hist->Integral();
      double Dybbar= Dybbar_hist->Integral();

      cout<<"In "+hist_name+", Data = "<<Data<<", All MC = "<<allMC<<", DY+b = "<<Dyb<<", DY+bbar = "<<Dybbar<<", and Charge = "<<(charges.at(ch)=="P"? Dybbar: Dyb)/(Dyb+Dybbar)<<", and Purity = "<<(Dyb+Dybbar)/allMC<<", and Norm = "<<Data/allMC<<endl;
    }
  }
}

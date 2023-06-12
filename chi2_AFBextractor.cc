#include "canvas_margin.h"
#include"AFBPlotter.cc"

AFBPlotter a;
TH1* AFB;
vector<double> sth2_values = {0.23151, 0.23154, 0.23157, 0.2230, 0.2300, 0.2305, 0.2310, 0.2315, 0.2320, 0.2325, 0.2330 };
void Hists_1D_AFB(TString channel="[em][em]bx201[6-8][ab]?/", TString charge="old", TString frame="AFBrecoil");
void Hists_2D_AFB(TString channel="[em][em]bx201[6-8][ab]?/", TString frame="AFBrecoil");
void Plots_1D_chi2(TString inputfile="sintheta.root");
void Plots_2D_chi2(TString inputfile="sintheta.root");


void chi2_AFBextractor(){
  //Hists_1D_AFB();
  //Hists_2D_AFB();

  //Plots_1D_chi2();
  Plots_2D_chi2();
}

void Hists_1D_AFB(TString channel="[em][em]bx201[6-8][ab]?/", TString charge="old", TString frame="AFBrecoil"){
  TString bjetcharge = "_[PM][0-5](x)";
  if(charge == "old") bjetcharge = "_[PM][2-5](x)";

  AFB = a.GetHist(0, channel+frame+bjetcharge, "");
  AFB->SetName("Data_sthw2");
  AFB->SaveAs("sthw2_Data.root");
  AFB = a.GetHist(1, channel+frame+bjetcharge, "");
  AFB->SetName("MC_sthw2");
  AFB->SaveAs("sthw2_MC.root");

  for(unsigned int i=0; i<sth2_values.size(); i++){
    AFB = a.GetHist(1, channel+frame+bjetcharge, Form("suffix:_sthw2_%d:dy",i));
    AFB->SetName(Form("MC_sthw2_%d",i));
    AFB->SaveAs(Form("sthw2_%d_MC.root",i));
  }
}

void Hists_2D_AFB(TString channel="[em][em]bx201[6-8][ab]?/", TString frame="AFBrecoil"){
  for(unsigned int ch=0; ch<6; ch++){
    AFB = a.GetHist(0, Form(channel+frame+"_[PM]%d(x)",ch), "");
    AFB->SetName(Form("Data_sthw2_ch%d",ch));
    AFB->SaveAs(Form("sthw2_ch%d_Data.root",ch));
    AFB = a.GetHist(1, Form(channel+frame+"_[PM]%d(x)",ch), "");
    AFB->SetName(Form("MC_sthw2_ch%d",ch));
    AFB->SaveAs(Form("sthw2_ch%d_MC.root",ch));

    for(unsigned int i=0; i<sth2_values.size(); i++){
      AFB = a.GetHist(1, Form(channel+frame+"_[PM]%d(x)",ch), Form("suffix:_sthw2_%d:dy",i));
      AFB->SetName(Form("MC_sthw2_%d_ch%d",i,ch));
      AFB->SaveAs(Form("sthw2_%d_ch%d_MC.root",i,ch));
    }
  }
}

void Plots_1D_chi2(TString inputfile="sintheta.root"){

  TFile *file_AFB = new TFile(inputfile);

  TCanvas *c_AFB = new TCanvas("c_AFB", "", 1000, 1000);
  c_AFB->Draw();

  TH1D *hist_data = (TH1D*)file_AFB->Get("Data_sthw2");
  TH1D *hist_mc   = (TH1D*)file_AFB->Get("MC_sthw2");

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
    TH1D *hist_mc_sthw2 = (TH1D*)file_AFB->Get(Form("MC_sthw2_%d",i));
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
  c_AFB->SaveAs("AFBrecoil_plots.png");
}

void Plots_2D_chi2(TString inputfile="sintheta.root"){

  TFile *file_AFB = new TFile(inputfile);

  for(unsigned int ch=0; ch<6; ch++){

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
    c_AFB->SaveAs(Form("AFBrecoil_plots_ch%d.png",ch));
  }
}

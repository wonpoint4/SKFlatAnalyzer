void Z0Weight(int mode = 0, TString Analyzer = "default"){

  AFBPlotter a;
  a.Reset();
  vector<TString> year = {"2016", "2017", "2018"};

  if(Analyzer == "default"){
    cout<<endl<<"Getting Z0Weight using preliminary results of GetZ0Weight in /data9/Users/wonjun/public/GetZ0Weight/"<<endl<<endl;;
    a.ScanFiles("/data9/Users/wonjun/public/GetZ0Weight/");
    Analyzer = "GetZ0Weight";
  }
  else{
    a.ScanFiles(TString()+getenv("SKFlatOutputDir")+getenv("SKFlatV")+"/"+Analyzer+"/");
  }
  
  samples["data"]=Sample("data",Sample::Type::DATA,kBlack,20)+TRegexp("201[0-9]/DATA/"+Analyzer+"_SkimTree_Dilepton_DoubleMuon_[A-Z]")+TRegexp("201[0-9]/DATA/"+Analyzer+"_SkimTree_Dilepton_DoubleEG_[A-Z]")+TRegexp("201[0-9]/DATA/"+Analyzer+"_SkimTree_Dilepton_MuonEG_[A-Z]")+TRegexp("201[0-9]/DATA/"+Analyzer+"_SkimTree_Dilepton_EGamma_[A-Z]");
  samples["mm2016"]=Sample("data (#mu#mu2016)",Sample::Type::DATA,kBlack,20)+TRegexp("2016/DATA/"+Analyzer+"_SkimTree_Dilepton_DoubleMuon_[A-Z]");
  samples["mm2017"]=Sample("data (#mu#mu2017)",Sample::Type::DATA,kRed,20)+TRegexp("2017/DATA/"+Analyzer+"_SkimTree_Dilepton_DoubleMuon_[A-Z]");
  samples["mm2018"]=Sample("data (#mu#mu2018)",Sample::Type::DATA,kBlue,20)+TRegexp("2018/DATA/"+Analyzer+"_SkimTree_Dilepton_DoubleMuon_[A-Z]");
  samples["ee2016"]=Sample("data (ee2016)",Sample::Type::DATA,kBlack,22)+TRegexp("2016/DATA/"+Analyzer+"_SkimTree_Dilepton_DoubleEG_[A-Z]");
  samples["ee2017"]=Sample("data (ee2017)",Sample::Type::DATA,kRed,22)+TRegexp("2017/DATA/"+Analyzer+"_SkimTree_Dilepton_DoubleEG_[A-Z]");
  samples["ee2018"]=Sample("data (ee2018)",Sample::Type::DATA,kBlue,22)+TRegexp("2018/DATA/"+Analyzer+"_SkimTree_Dilepton_EGamma_[A-Z]");
  samples["em2016"]=Sample("data (e#mu2016)",Sample::Type::DATA,kBlack,22)+TRegexp("2016/DATA/"+Analyzer+"_SkimTree_Dilepton_MuonEG_[A-Z]");
  samples["em2017"]=Sample("data (e#mu2017)",Sample::Type::DATA,kRed,22)+TRegexp("2017/DATA/"+Analyzer+"_SkimTree_Dilepton_MuonEG_[A-Z]");
  samples["em2018"]=Sample("data (e#mu2018)",Sample::Type::DATA,kBlue,22)+TRegexp("2018/DATA/"+Analyzer+"_SkimTree_Dilepton_MuonEG_[A-Z]");
  samples["mm"]=Sample("data (#mu#mu)",Sample::Type::DATA,kBlack,20)+TRegexp("201[0-9]/DATA/"+Analyzer+"_SkimTree_Dilepton_DoubleMuon_[A-Z]");
  samples["ee"]=Sample("data (ee)",Sample::Type::DATA,kRed,22)+TRegexp("201[0-9]/DATA/"+Analyzer+"_SkimTree_Dilepton_DoubleEG_[A-Z]")+TRegexp("201[0-9]/DATA/"+Analyzer+"_SkimTree_Dilepton_EGamma_[A-Z]");
  samples["em"]=Sample("data (e#mu)",Sample::Type::DATA,kRed,22)+TRegexp("201[0-9]/DATA/"+Analyzer+"_SkimTree_Dilepton_MuonEG_[A-Z]");

  samples["mc"]=Sample("mc",Sample::Type::SIGNAL,kBlue)+TRegexp("/"+Analyzer+"_SkimTree_Dilepton_DYJets$")+TRegexp("/"+Analyzer+"_SkimTree_Dilepton_[W-Z][W-Z]_pythia")+TRegexp("/"+Analyzer+"_SkimTree_Dilepton_WJets_MG")+TRegexp("/"+Analyzer+"_SkimTree_Dilepton_TTLL_powheg")+TRegexp("/*/"+Analyzer+"_TTLJ_powheg");
  samples["ww"]=Sample("WW",Sample::Type::BG,kBlue)+TRegexp("/"+Analyzer+"_SkimTree_Dilepton_WW_pythia");
  samples["wz"]=Sample("WZ",Sample::Type::BG,kGreen)+TRegexp("/"+Analyzer+"_SkimTree_Dilepton_WZ_pythia");
  samples["zz"]=Sample("ZZ",Sample::Type::BG,kCyan)+TRegexp("/"+Analyzer+"_SkimTree_Dilepton_ZZ_pythia");
  samples["vv"]=Sample("Diboson",Sample::Type::BG,kBlue)+TRegexp("/"+Analyzer+"_SkimTree_Dilepton_[W-Z][W-Z]_pythia");
  samples["wjets"]=Sample("W",Sample::Type::BG,kYellow)+TRegexp("/"+Analyzer+"_SkimTree_Dilepton_WJets_MG");
  samples["tt"]=Sample("t#bar{t}",Sample::Type::BG,kMagenta)+TRegexp("/"+Analyzer+"_SkimTree_Dilepton_TTLL_powheg");
  samples["ttlj"]=Sample("TTLJ",Sample::Type::BG,kMagenta+6)+TRegexp("/*/"+Analyzer+"_TTLJ_powheg");
  samples["amc"]=Sample("#gamma*/Z#rightarrowll",Sample::Type::SIGNAL,kRed)+TRegexp("/"+Analyzer+"_SkimTree_Dilepton_DYJets$");  
  samples["tau_amc"]="tau_"%(Sample("#gamma*/Z#rightarrow#tau#tau",Sample::Type::BG,kGreen)+"amc");

  //mode == 0 :Getting Z0Weight
  if(mode == 0 ){
    a.SetupEntries("data mc");
    TString command = "";
    for(unsigned int i=0;i<year.size();i++){
      TFile file("Z0Weight_won_"+year.at(i)+".root","recreate");

      TH1* Data = a.GetHist("data", "ll"+year.at(i)+"/z0_noz0");
      if(!Data){
	cout<<"Empty histograms of data, ll"+year.at(i)+"/z0_noz0"<<endl;
	continue;
      }
      TH1* MC = a.GetHist("mc", "ll"+year.at(i)+"/z0_noz0");
      if(!MC){
	cout<<"Empty histograms of mc, ll"+year.at(i)+"/z0_noz0"<<endl;
	continue;
      }
      TF1* Data_fit = new TF1("data_fit", "gaus", -10, 10);
      TF1* MC_fit = new TF1("mc_fit","gaus",-10, 10);

      Data->SetTitle("ll"+year.at(i)+"/z0[-10,10]/z0_noz0_Data");
      MC->SetTitle("ll"+year.at(i)+"/z0[-10,10]/z0_noz0_MC");
      Data_fit->SetTitle("ll"+year.at(i)+"/z0[-10,10]/z0_noz0_Data_GaussianFit");
      MC_fit->SetTitle("ll"+year.at(i)+"/z0[-10,10]/z0_noz0_MC_GaussianFit");

      Data->Fit(Data_fit,"","",-10,10);
      MC->Fit(MC_fit,"","",-10,10);

      //In case of getting z0weight from another phase space Ex) y[-2,2] or z0[-15,15]

      //TH1* Data_central = a.GetHist("data", "ll"+year.at(i)+"/z0_noz0");
      //TH1* MC_central = a.GetHist("mc", "ll"+year.at(i)+"/z0_noz0");
      //TF1* Data_fit_central = new TF1("data_fit_z0_less15", "gaus", -15, 15);
      //TF1* MC_fit_central = new TF1("mc_fit_z0_less15","gaus",-15, 15);
      //Data_central->GetXaxis()->SetRangeUser(-15,15);
      //MC_central->GetXaxis()->SetRangeUser(-15,15);

      //Data_central->SetTitle("ll"+year.at(i)+"/z0[-15,15]/z0_noz0_Data");
      //MC_central->SetTitle("ll"+year.at(i)+"/z0[-15,15]/z0_noz0_MC");
      //Data_fit_central->SetTitle("ll"+year.at(i)+"/z0[-15,15]/z0_noz0_Data");
      //MC_fit_central->SetTitle("ll"+year.at(i)+"/z0[-15,15]/z0_noz0_MC");

      //Data_central->Fit(Data_fit_central);
      //MC_central->Fit(MC_fit_central);

      file.cd();
      Data->Write();
      MC->Write();
      Data_fit->Write();
      MC_fit->Write();
      //Data_fit_central->Write();
      //MC_fit_central->Write();

      file.Close();
      command += "mv Z0Weight_won_"+year.at(i)+".root data/Run2Legacy_v4/"+year.at(i)+"/Z0/Z0Weight.root\n";
    }
    cout<<endl<<command<<endl;
  }

  //mode == 1 : Drawing plots for AN
  else if(mode == 1){
    for(unsigned int i=0;i<year.size();i++){
      vector<TString> channel = {"mm","ee", "ll"};
      for(unsigned int j=0;j<channel.size();j++){
	a.SetupEntries("data ^amc+tau_amc+vv+wjets+tt+ttlj");
	a.SavePlot(channel.at(j)+year.at(i)+"_z0_beforeZ0",        "histname:"+channel.at(j)+year.at(i)+"/z0_noz0 norm pdf");
	a.SavePlot(channel.at(j)+year.at(i)+"_z0_afterZ0",         "histname:"+channel.at(j)+year.at(i)+"/z0 norm pdf");
	a.SavePlot(channel.at(j)+year.at(i)+"_l1eta_beforeZ0",     "histname:"+channel.at(j)+year.at(i)+"/l1eta_noz0 norm pdf");
        a.SavePlot(channel.at(j)+year.at(i)+"_l1eta_afterZ0",      "histname:"+channel.at(j)+year.at(i)+"/l1eta norm pdf");
	//a.SavePlot(channel.at(j)+year.at(i)+"_dirap_beforeZ0",     "histname:"+channel.at(j)+year.at(i)+"/dirap_noz0 norm pdf"); //dirap seems doesn't work in SavePlot
        //a.SavePlot(channel.at(j)+year.at(i)+"_dirap_afterZ0",      "histname:"+channel.at(j)+year.at(i)+"/dirap norm pdf");
	a.DrawPlot(channel.at(j)+year.at(i)+"/dirap_noz0", "norm")->SaveAs(channel.at(j)+year.at(i)+"_dirap_beforeZ0.pdf");
	a.DrawPlot(channel.at(j)+year.at(i)+"/dirap", "norm")->SaveAs(channel.at(j)+year.at(i)+"_dirap_afterZ0.pdf");
        a.SavePlot(channel.at(j)+year.at(i)+"_AFB(dirap)_beforeZ0","histname:"+channel.at(j)+year.at(i)+"/AFB_noz0(y) norm pdf");
        a.SavePlot(channel.at(j)+year.at(i)+"_AFB(dirap)_afterZ0", "histname:"+channel.at(j)+year.at(i)+"/AFB(y) norm pdf");
      }
    }
  }
  
  //modes : Comparison normalization of MC before/after z0 weight
  else{
    a.SetupEntries("data ^amc+tau_amc+vv+wjets+tt+ttlj");
    for(unsigned int i=0;i<year.size();i++){
      vector<TString> corrected = {"_", "_nocos_", "_nozpt_", "_noSF_", "_noPU_"};
      vector<TString> corrected2 = {"", "_nocos", "_nozpt", "_noSF", "_noPU"};

      for(unsigned int j=0; j<corrected.size(); j++){
	long int1 = a.GetHist("mc", "ll"+year.at(i)+"/z0"+corrected.at(j)+"noz0")->Integral(1,120);
	long int2 = a.GetHist("mc", "ll"+year.at(i)+"/z0"+corrected.at(j)+"noz0")->Integral();
	long int3 = a.GetHist("mc", "ll"+year.at(i)+"/z0"+corrected2.at(j))->Integral(1,120);
	long int4 = a.GetHist("mc", "ll"+year.at(i)+"/z0"+corrected2.at(j))->Integral();
	//long int5 = a.GetHist("mc", "ll"+year.at(i)+"/z0"+corrected.at(j)+"fitz0")->Integral(1,120);
	//long int6 = a.GetHist("mc", "ll"+year.at(i)+"/z0"+corrected.at(j)+"fitz0")->Integral();
	//long int7 = a.GetHist("mc", "ll"+year.at(i)+"/z0"+corrected.at(j)+"fitz1")->Integral(1,120);
	//long int8 = a.GetHist("mc", "ll"+year.at(i)+"/z0"+corrected.at(j)+"fitz1")->Integral();

	cout<<"ll"+year.at(i)+"/z0"+corrected.at(j)+"noz0 integral(1,120) = "<<int1<<endl;
	cout<<"ll"+year.at(i)+"/z0"+corrected.at(j)+"noz0 integral(0,121) = "<<int2<<endl;
	cout<<"ll"+year.at(i)+"/z0"+corrected2.at(j)+" integral(1,120) = "<<int3<<", norm difference = "<<int3-int1<<endl;
	cout<<"ll"+year.at(i)+"/z0"+corrected2.at(j)+" integral(0,121) = "<<int4<<", norm difference = "<<int4-int2<<endl;
	//cout<<"ll"+year.at(i)+"/z0"+corrected.at(j)+"fitz0 integral(1,120) = "<<int5<<", norm difference = "<<int5-int1<<endl;
	//cout<<"ll"+year.at(i)+"/z0"+corrected.at(j)+"fitz0 integral(0,121) = "<<int6<<", norm difference = "<<int6-int2<<endl;
	//cout<<"ll"+year.at(i)+"/z0"+corrected.at(j)+"fitz1 integral(1,120) = "<<int7<<", norm difference = "<<int7-int1<<endl;
	//cout<<"ll"+year.at(i)+"/z0"+corrected.at(j)+"fitz1 integral(0,121) = "<<int8<<", norm difference = "<<int8-int2<<endl;
      }
    }
  }
}

{
  EfficiencyPlotter amc("data ^amc+tau_amc+vv+wjets+tttw");
  vector<TString> eras={"ee2016a","ee2016b","ee2017","ee2018"};
  Verbosity=0;
  for(auto era:eras){
    amc.entries[1].weight=1.;
    TH1* hdata=amc.GetHist(0,era+"/m60to120/dimass");
    TH1* hmc=amc.GetTH1(amc.GetHist(1,era+"/m60to120/dimass"));
    amc.entries[1].weight=hdata->Integral()/hmc->Integral();
    amc.SavePlot("210326/"+era+"_dimass","histname:"+era+"/m60to120/dimass widthweight 2:widewidey");
    amc.SavePlot("210326/"+era+"_ss_dimass","histname:"+era+"/m60to120/ss_dimass widthweight 2:widewidey");
    TH1* hdata_os=amc.GetHist(0,era+"/m60to120/dimass");
    TH1* hmc_os=amc.GetTH1(amc.GetHist(1,era+"/m60to120/dimass"));
    TH1* hdata_ss=amc.GetHist(0,era+"/m60to120/ss_dimass");
    TH1* hmc_ss=amc.GetTH1(amc.GetHist(1,era+"/m60to120/ss_dimass"));
    cout<<"Norm factor: "<<amc.entries[1].weight<<endl;
    cout<<"OS Integral: "<<hdata_os->Integral()<<" "<<hmc_os->Integral()<<endl;
    cout<<"SS Integral: "<<hdata_ss->Integral()<<" "<<hmc_ss->Integral()<<endl;
    cout<<"CF Rate: "<<hdata_ss->Integral()/hdata_os->Integral()/2<<" "<<hmc_ss->Integral()/hmc_os->Integral()/2<<endl;
    cout<<"CF SF: "<<hdata_ss->Integral()/hdata_os->Integral()/2/(hmc_ss->Integral()/hmc_os->Integral()/2)<<endl;
  }
  AFBPlotter aamc;
  aamc.SavePlot("210326/EE2017_dimass","histname:EE2017/m[60,3000]/dimass widthweight logx type:1");
  aamc.SavePlot("210326/EE2017_ss_dimass","histname:EE2017/m[60,3000]/ss_dimass widthweight logx type:1");
  aamc.SavePlot("210326/eE2017_dimass","histname:eE2017/m[60,3000]/dimass widthweight logx type:1");
  aamc.SavePlot("210326/eE2017_ss_dimass","histname:eE2017/m[60,3000]/ss_dimass widthweight logx type:1");

  aamc.SavePlot("210326/MM2017_dimass","histname:MM2017/m[60,3000]/dimass widthweight logx type:1");
  aamc.SavePlot("210326/MM2017_ss_dimass","histname:MM2017/m[60,3000]/ss_dimass widthweight logx type:1");
  aamc.SavePlot("210326/mM2017_dimass","histname:mM2017/m[60,3000]/dimass widthweight logx type:1");
  aamc.SavePlot("210326/mM2017_ss_dimass","histname:mM2017/m[60,3000]/ss_dimass widthweight logx type:1");

  aamc.Setup("data-amc-tau_amc-vv-wjets-tttw data-amc-tau_amc-vv-wjets-tttw");
  aamc.entries[1].SetStyle(2);
  aamc.entries[0].title="OS data-MC";
  aamc.entries[1].title="SS data-MC";
  aamc.entries[1].SetHistPrefix("ss_");
  aamc.SavePlot("210326/EE2017_dimass_sss","histname:EE2017/m[60,3000]/dimass widthweight logx 2:widewidey");
  aamc.SavePlot("210326/eE2017_dimass_sss","histname:eE2017/m[60,3000]/dimass widthweight logx 2:widewidey");
  aamc.SavePlot("210326/MM2017_dimass_sss","histname:MM2017/m[60,3000]/dimass widthweight logx 2:widewidey");
  aamc.SavePlot("210326/mM2017_dimass_sss","histname:mM2017/m[60,3000]/dimass widthweight logx 2:widewidey");

}
    

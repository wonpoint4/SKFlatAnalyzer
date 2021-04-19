void plot1(){
  EfficiencyPlotter amc("data ^amc+tau_amc+vv+wjets+tttw");
  EfficiencyPlotter amcss("data ^amc+tau_amc+vv+wjets+tttw+ss");
  EfficiencyPlotter amcss2("data ^amc+tau_amc+vv+wjets+tttw+1.7*ss");
  EfficiencyPlotter amcS20("data ^amcS20+tau_amcS20+vv+wjets+tttw");
  EfficiencyPlotter amcS20ss("data ^amcS20+tau_amcS20+vv+wjets+tttw+ss");
  EfficiencyPlotter amcS20ss2("data ^amcS20+tau_amcS20+vv+wjets+tttw+1.7*ss");
  EfficiencyPlotter mg("data ^mg+tau_mg+vv+wjets+tttw");
  EfficiencyPlotter mgss("data ^mg+tau_mg+vv+wjets+tttw+ss");
  EfficiencyPlotter mgss2("data ^mg+tau_mg+vv+wjets+tttw+1.7*ss");
  EfficiencyPlotter mi("data ^minnlo+tau_minnlo+vv+wjets+tttw");
  EfficiencyPlotter miss("data ^minnlo+tau_minnlo+vv+wjets+tttw+ss");
  EfficiencyPlotter miss2("data ^minnlo+tau_minnlo+vv+wjets+tttw+1.7*ss");
  for(TString channel:{"mm","ee"}){
    for(TString era:{"2016a","2016b","2017","2018"}){
      TString title="#mu#mu";
      if(channel=="ee") title="ee";
      if(era=="2016a") title+=" 2016preVFP";
      else if(era=="2016b") title+=" 2016postVFP";
      else title+=" "+era;
      TString lepton="Muon";
      if(channel=="ee") lepton="Electron";
      TString ll="#mu#mu";
      if(channel=="ee") ll="ee";
      amc.SavePlot("KPS/eff/"+channel+era+"_dimass_amc","histname:"+channel+era+"/m52to150/dimass norm rebin:2 1:logy xtitle:'m("+ll+") [GeV]' title++ title:'"+title+"'");
      amc.SavePlot("KPS/eff/"+channel+era+"_leta_amc","histname:"+channel+era+"/m52to150/leta norm 1:logy xmin:-2.5 xmax:2.5 noleg xtitle:'"+lepton+" #eta' title++ title:'"+title+"'");
      amc.SavePlot("KPS/eff/"+channel+era+"_lpt_amc","histname:"+channel+era+"/m52to150/lpt norm 1:logy xmax:100 xtitle:'"+lepton+" p_{T} [GeV]' title++ title:'"+title+"'");
      if(channel.Contains("e")){
	amcss.SavePlot("KPS/eff/"+channel+era+"_dimass_amcss","histname:"+channel+era+"/m52to150/dimass norm rebin:2 1:logy xtitle:'m("+ll+") [GeV]' title++ title:'"+title+"'");
	amcss.SavePlot("KPS/eff/"+channel+era+"_leta_amcss","histname:"+channel+era+"/m52to150/leta norm 1:logy xmin:-2.5 xmax:2.5 noleg xtitle:'"+lepton+" #eta' title++ title:'"+title+"'");
	amcss.SavePlot("KPS/eff/"+channel+era+"_lpt_amcss","histname:"+channel+era+"/m52to150/lpt norm 1:logy xmax:100 xtitle:'"+lepton+" p_{T} [GeV]' title++ title:'"+title+"'");
	amcss.SavePlot("KPS/eff/"+channel+era+"_lsceta_amcss","histname:"+channel+era+"/m52to150/lsceta norm 1:logy xmin:-2.5 xmax:2.5 noleg xtitle:'"+lepton+" #eta_{SC}' title++ title:'"+title+"'");
	amc.SavePlot("KPS/eff/"+channel+era+"_lsceta_amc","histname:"+channel+era+"/m52to150/lsceta norm 1:logy xmin:-2.5 xmax:2.5 noleg xtitle:'"+lepton+" #eta_{SC}' title++ title:'"+title+"'");
      }
      if(channel.Contains("m")){
	amcss2.SavePlot("KPS/eff/"+channel+era+"_dimass_amcss2","histname:"+channel+era+"/m52to150/dimass norm rebin:2 1:logy xtitle:'m("+ll+") [GeV]' title++ title:'"+title+"'");
	amcss2.SavePlot("KPS/eff/"+channel+era+"_leta_amcss2","histname:"+channel+era+"/m52to150/leta norm 1:logy xmin:-2.5 xmax:2.5 noleg xtitle:'"+lepton+" #eta' title++ title:'"+title+"'");
	amcss2.SavePlot("KPS/eff/"+channel+era+"_lpt_amcss2","histname:"+channel+era+"/m52to150/lpt norm 1:logy xmax:100 xtitle:'"+lepton+" p_{T} [GeV]' title++ title:'"+title+"'");
	amcss2.SavePlot("KPS/eff/"+channel+era+"_lpt_zpeak_amcss2","histname:"+channel+era+"/m80to100/lpt norm 1:logy xmax:100 xtitle:'"+lepton+" p_{T} [GeV]' title++ title:'"+title+"'");
      }

      mg.SavePlot("KPS/eff/"+channel+era+"_dimass_mg","histname:"+channel+era+"/m52to150/dimass norm rebin:2 1:logy xtitle:'m("+ll+") [GeV]' title++ title:'"+title+"'");
      mg.SavePlot("KPS/eff/"+channel+era+"_leta_mg","histname:"+channel+era+"/m52to150/leta norm 1:logy xmin:-2.5 xmax:2.5 noleg xtitle:'"+lepton+" #eta' title++ title:'"+title+"'");
      mg.SavePlot("KPS/eff/"+channel+era+"_lpt_mg","histname:"+channel+era+"/m52to150/lpt norm 1:logy xmax:100 xtitle:'"+lepton+" p_{T} [GeV]' title++ title:'"+title+"'");
      if(channel.Contains("e")){
	mgss.SavePlot("KPS/eff/"+channel+era+"_dimass_mgss","histname:"+channel+era+"/m52to150/dimass norm rebin:2 1:logy xtitle:'m("+ll+") [GeV]' title++ title:'"+title+"'");
	mgss.SavePlot("KPS/eff/"+channel+era+"_leta_mgss","histname:"+channel+era+"/m52to150/leta norm 1:logy xmin:-2.5 xmax:2.5 noleg xtitle:'"+lepton+" #eta' title++ title:'"+title+"'");
	mgss.SavePlot("KPS/eff/"+channel+era+"_lpt_mgss","histname:"+channel+era+"/m52to150/lpt norm 1:logy xmax:100 xtitle:'"+lepton+" p_{T} [GeV]' title++ title:'"+title+"'");
	mgss.SavePlot("KPS/eff/"+channel+era+"_lsceta_mgss","histname:"+channel+era+"/m52to150/lsceta norm 1:logy xmin:-2.5 xmax:2.5 noleg xtitle:'"+lepton+" #eta_{SC}' title++ title:'"+title+"'");
	mg.SavePlot("KPS/eff/"+channel+era+"_lsceta_mg","histname:"+channel+era+"/m52to150/lsceta norm 1:logy xmin:-2.5 xmax:2.5 noleg xtitle:'"+lepton+" #eta_{SC}' title++ title:'"+title+"'");
      }
      if(channel.Contains("m")){
	mgss2.SavePlot("KPS/eff/"+channel+era+"_dimass_mgss2","histname:"+channel+era+"/m52to150/dimass norm rebin:2 1:logy xtitle:'m("+ll+") [GeV]' title++ title:'"+title+"'");
	mgss2.SavePlot("KPS/eff/"+channel+era+"_leta_mgss2","histname:"+channel+era+"/m52to150/leta norm 1:logy xmin:-2.5 xmax:2.5 noleg xtitle:'"+lepton+" #eta' title++ title:'"+title+"'");
	mgss2.SavePlot("KPS/eff/"+channel+era+"_lpt_mgss2","histname:"+channel+era+"/m52to150/lpt norm 1:logy xmax:100 xtitle:'"+lepton+" p_{T} [GeV]' title++ title:'"+title+"'");
	mgss2.SavePlot("KPS/eff/"+channel+era+"_lpt_zpeak_mgss2","histname:"+channel+era+"/m80to100/lpt norm 1:logy xmax:100 xtitle:'"+lepton+" p_{T} [GeV]' title++ title:'"+title+"'");
      }

      if(!era.Contains("2016")) continue;
      mi.SavePlot("KPS/eff/"+channel+era+"_dimass_mi","histname:"+channel+era+"/m52to150/dimass norm rebin:2 1:logy xtitle:'m("+ll+") [GeV]' title++ title:'"+title+"'");
      mi.SavePlot("KPS/eff/"+channel+era+"_leta_mi","histname:"+channel+era+"/m52to150/leta norm 1:logy xmin:-2.5 xmax:2.5 noleg xtitle:'"+lepton+" #eta' title++ title:'"+title+"'");
      mi.SavePlot("KPS/eff/"+channel+era+"_lpt_mi","histname:"+channel+era+"/m52to150/lpt norm 1:logy xmax:100 xtitle:'"+lepton+" p_{T} [GeV]' title++ title:'"+title+"'");
      if(channel.Contains("e")){
	miss.SavePlot("KPS/eff/"+channel+era+"_dimass_miss","histname:"+channel+era+"/m52to150/dimass norm rebin:2 1:logy xtitle:'m("+ll+") [GeV]' title++ title:'"+title+"'");
	miss.SavePlot("KPS/eff/"+channel+era+"_leta_miss","histname:"+channel+era+"/m52to150/leta norm 1:logy xmin:-2.5 xmax:2.5 noleg xtitle:'"+lepton+" #eta' title++ title:'"+title+"'");
	miss.SavePlot("KPS/eff/"+channel+era+"_lpt_miss","histname:"+channel+era+"/m52to150/lpt norm 1:logy xmax:100 xtitle:'"+lepton+" p_{T} [GeV]' title++ title:'"+title+"'");
	miss.SavePlot("KPS/eff/"+channel+era+"_lsceta_miss","histname:"+channel+era+"/m52to150/lsceta norm 1:logy xmin:-2.5 xmax:2.5 noleg xtitle:'"+lepton+" #eta_{SC}' title++ title:'"+title+"'");
	mi.SavePlot("KPS/eff/"+channel+era+"_lsceta_mi","histname:"+channel+era+"/m52to150/lsceta norm 1:logy xmin:-2.5 xmax:2.5 noleg xtitle:'"+lepton+" #eta_{SC}' title++ title:'"+title+"'");
      }
      if(channel.Contains("m")){
	miss2.SavePlot("KPS/eff/"+channel+era+"_dimass_miss2","histname:"+channel+era+"/m52to150/dimass norm rebin:2 1:logy xtitle:'m("+ll+") [GeV]' title++ title:'"+title+"'");
	miss2.SavePlot("KPS/eff/"+channel+era+"_leta_miss2","histname:"+channel+era+"/m52to150/leta norm 1:logy xmin:-2.5 xmax:2.5 noleg xtitle:'"+lepton+" #eta' title++ title:'"+title+"'");
	miss2.SavePlot("KPS/eff/"+channel+era+"_lpt_miss2","histname:"+channel+era+"/m52to150/lpt norm 1:logy xmax:100 xtitle:'"+lepton+" p_{T} [GeV]' title++ title:'"+title+"'");
	miss2.SavePlot("KPS/eff/"+channel+era+"_lpt_zpeak_miss2","histname:"+channel+era+"/m80to100/lpt norm 1:logy xmax:100 xtitle:'"+lepton+" p_{T} [GeV]' title++ title:'"+title+"'");
      }
      if(!era.Contains("2016a")) continue;
      amcS20.SavePlot("KPS/eff/"+channel+era+"_dimass_amcS20","histname:"+channel+era+"/m52to150/dimass norm rebin:2 1:logy xtitle:'m("+ll+") [GeV]' title++ title:'"+title+"'");
      amcS20.SavePlot("KPS/eff/"+channel+era+"_leta_amcS20","histname:"+channel+era+"/m52to150/leta norm 1:logy xmin:-2.5 xmax:2.5 noleg xtitle:'"+lepton+" #eta' title++ title:'"+title+"'");
      amcS20.SavePlot("KPS/eff/"+channel+era+"_lpt_amcS20","histname:"+channel+era+"/m52to150/lpt norm 1:logy xmax:100 xtitle:'"+lepton+" p_{T} [GeV]' title++ title:'"+title+"'");
      if(channel.Contains("e")){
	amcS20ss.SavePlot("KPS/eff/"+channel+era+"_dimass_amcS20ss","histname:"+channel+era+"/m52to150/dimass norm rebin:2 1:logy xtitle:'m("+ll+") [GeV]' title++ title:'"+title+"'");
	amcS20ss.SavePlot("KPS/eff/"+channel+era+"_leta_amcS20ss","histname:"+channel+era+"/m52to150/leta norm 1:logy xmin:-2.5 xmax:2.5 noleg xtitle:'"+lepton+" #eta' title++ title:'"+title+"'");
	amcS20ss.SavePlot("KPS/eff/"+channel+era+"_lpt_amcS20ss","histname:"+channel+era+"/m52to150/lpt norm 1:logy xmax:100 xtitle:'"+lepton+" p_{T} [GeV]' title++ title:'"+title+"'");
	amcS20ss.SavePlot("KPS/eff/"+channel+era+"_lsceta_amcS20ss","histname:"+channel+era+"/m52to150/lsceta norm 1:logy xmin:-2.5 xmax:2.5 noleg xtitle:'"+lepton+" #eta_{SC}' title++ title:'"+title+"'");
	amcS20.SavePlot("KPS/eff/"+channel+era+"_lsceta_amcS20","histname:"+channel+era+"/m52to150/lsceta norm 1:logy xmin:-2.5 xmax:2.5 noleg xtitle:'"+lepton+" #eta_{SC}' title++ title:'"+title+"'");
      }
      if(channel.Contains("m")){
	amcS20ss2.SavePlot("KPS/eff/"+channel+era+"_dimass_amcS20ss2","histname:"+channel+era+"/m52to150/dimass norm rebin:2 1:logy xtitle:'m("+ll+") [GeV]' title++ title:'"+title+"'");
	amcS20ss2.SavePlot("KPS/eff/"+channel+era+"_leta_amcS20ss2","histname:"+channel+era+"/m52to150/leta norm 1:logy xmin:-2.5 xmax:2.5 noleg xtitle:'"+lepton+" #eta' title++ title:'"+title+"'");
	amcS20ss2.SavePlot("KPS/eff/"+channel+era+"_lpt_amcS20ss2","histname:"+channel+era+"/m52to150/lpt norm 1:logy xmax:100 xtitle:'"+lepton+" p_{T} [GeV]' title++ title:'"+title+"'");
	amcS20ss2.SavePlot("KPS/eff/"+channel+era+"_lpt_zpeak_amcS20ss2","histname:"+channel+era+"/m80to100/lpt norm 1:logy xmax:100 xtitle:'"+lepton+" p_{T} [GeV]' title++ title:'"+title+"'");
      }
    }
  }
}
void plot2(){
  EfficiencyPlotter amc("data ^amc+tau_amc+vv+wjets+tttw");
  amc.SavePlot("KPS/eff/ee2017_l0pt_amc","histname:ee2017/m52to150/l0pt norm 1:logy xmax:100 xtitle:'Leading electron p_{T} [GeV]' title++ title:'ee 2017'");
  amc.SavePlot("KPS/eff/ee2017_l0pt_amc_noL1","histname:ee2017/m52to150/l0pt_noL1 norm 1:logy xmax:100 xtitle:'Leading electron p_{T} [GeV]' title++ title:'ee 2017'");
  amc.SavePlot("KPS/eff/ee2018_l0pt_amc","histname:ee2018/m52to150/l0pt norm 1:logy xmax:100 xtitle:'Leading electron p_{T} [GeV]' title++ title:'ee 2018'");
  amc.SavePlot("KPS/eff/ee2018_l0pt_amc_noL1","histname:ee2018/m52to150/l0pt_noL1 norm 1:logy xmax:100 xtitle:'Leading electron p_{T} [GeV]' title++ title:'ee 2018'");
}
void plot3(){
  AFBPlotter amcss("amc+tau_amc+vv+wjets+tttw ss");
  amcss.entries[1].style.linewidth=2;
  AFBPlotter amcss2("amc+tau_amc+vv+wjets+tttw 1.7*ss");
  amcss2.entries[1].style.linewidth=2;
  for(TString channel:{"mm","ee"}){
    for(TString era:{"2016a","2016b","2017","2018"}){
      TString title="#mu#mu";
      if(channel=="ee") title="ee";
      if(era=="2016a") title+=" 2016preVFP";
      else if(era=="2016b") title+=" 2016postVFP";
      else title+=" "+era;
      if(channel=="mm"){
	amcss2.SavePlot("KPS/qcd/"+channel+era+"_ratio","histname:"+channel+era+"/costhetaCS project:x xmin:52 xmax:1000 ymin:0 ymax:0.1 base:0 ytitle:'QCD multi-jet fraction' noleg logx type:2 title++ title:'"+title+"'");
      }else{
	amcss.SavePlot("KPS/qcd/"+channel+era+"_ratio","histname:"+channel+era+"/costhetaCS project:x xmin:52 xmax:1000 ymin:0 ymax:0.1 base:0 ytitle:'QCD multi-jet fraction' noleg logx type:2 title++ title:'"+title+"'");
      }
    }
  }
}
void plot4(){
  TCanvas* c;
  TFile* f;
  TH1 *hdata,*hmc,*axisowner;
  TLegend *leg;

  f=new TFile("data/Run2UltraLegacy_v1/2018/SMP/CFRate_old0.root");
  c=(TCanvas*)f->Get("cf_pt");
  leg=new TLegend(0.17,0.88,0.5,0.6);
  axisowner=(TH1*)c->GetPrimitive("cfdata_py_0.0_1.0");
  axisowner->GetXaxis()->SetTitle("Electron p_{T} [GeV]");
  axisowner->GetYaxis()->SetTitle("Charg-flip rate");
  ((TH1*)c->GetPrimitive("cfdata_py_2.0_2.5"))->SetLineColor(6);
  ((TH1*)c->GetPrimitive("cfmc_py_2.0_2.5"))->SetLineColor(6);
  hdata=new TH1D;
  hmc=new TH1D;
  hmc->SetLineStyle(2);
  leg->AddEntry(hdata,"data","l");
  leg->AddEntry(hmc,"MC","l");
  leg->AddEntry(c->GetPrimitive("cfdata_py_0.0_1.0"),"0.0 #leq |#eta| < 1.0","lp");
  leg->AddEntry(c->GetPrimitive("cfdata_py_1.0_1.4"),"1.0 #leq |#eta| < 1.4","lp");
  leg->AddEntry(c->GetPrimitive("cfdata_py_1.4_1.7"),"1.4 #leq |#eta| < 1.7","lp");
  leg->AddEntry(c->GetPrimitive("cfdata_py_1.7_2.0"),"1.7 #leq |#eta| < 2.0","lp");
  leg->AddEntry(c->GetPrimitive("cfdata_py_2.0_2.5"),"2.0 #leq |#eta| < 2.5","lp");
  leg->Draw();
  axisowner->SetTitle("2018");
  c->Draw();
  c->Update();
  ((TPaveText*)c->GetPrimitive("title"))->SetTextSize(0.075);
  c->SaveAs("fig/KPS/CFRate/ee2018.png");

  f=new TFile("data/Run2UltraLegacy_v1/2017/SMP/CFRate_old0.root");
  c=(TCanvas*)f->Get("cf_pt");
  leg=new TLegend(0.17,0.88,0.5,0.6);
  axisowner=(TH1*)c->GetPrimitive("cfdata_py_0.0_1.0");
  axisowner->GetXaxis()->SetTitle("Electron p_{T} [GeV]");
  axisowner->GetYaxis()->SetTitle("Charg-flip rate");
  axisowner->GetYaxis()->SetRangeUser(0,0.06);
  ((TH1*)c->GetPrimitive("cfdata_py_2.0_2.5"))->SetLineColor(6);
  ((TH1*)c->GetPrimitive("cfmc_py_2.0_2.5"))->SetLineColor(6);
  hdata=new TH1D;
  hmc=new TH1D;
  hmc->SetLineStyle(2);
  leg->AddEntry(hdata,"data","l");
  leg->AddEntry(hmc,"MC","l");
  leg->AddEntry(c->GetPrimitive("cfdata_py_0.0_1.0"),"0.0 #leq |#eta| < 1.0","lp");
  leg->AddEntry(c->GetPrimitive("cfdata_py_1.0_1.4"),"1.0 #leq |#eta| < 1.4","lp");
  leg->AddEntry(c->GetPrimitive("cfdata_py_1.4_1.7"),"1.4 #leq |#eta| < 1.7","lp");
  leg->AddEntry(c->GetPrimitive("cfdata_py_1.7_2.0"),"1.7 #leq |#eta| < 2.0","lp");
  leg->AddEntry(c->GetPrimitive("cfdata_py_2.0_2.5"),"2.0 #leq |#eta| < 2.5","lp");
  leg->Draw();
  axisowner->SetTitle("2017");
  c->Draw();
  c->Update();
  ((TPaveText*)c->GetPrimitive("title"))->SetTextSize(0.075);
  c->SaveAs("fig/KPS/CFRate/ee2017.png");
  
  f=new TFile("data/Run2UltraLegacy_v1/2016preVFP/SMP/CFRate_old0.root");
  c=(TCanvas*)f->Get("cf_pt");
  leg=new TLegend(0.17,0.88,0.5,0.6);
  axisowner=(TH1*)c->GetPrimitive("cfdata_py_0.0_1.0");
  axisowner->GetXaxis()->SetTitle("Electron p_{T} [GeV]");
  axisowner->GetYaxis()->SetTitle("Charg-flip rate");
  ((TH1*)c->GetPrimitive("cfdata_py_2.0_2.5"))->SetLineColor(6);
  ((TH1*)c->GetPrimitive("cfmc_py_2.0_2.5"))->SetLineColor(6);
  hdata=new TH1D;
  hmc=new TH1D;
  hmc->SetLineStyle(2);
  leg->AddEntry(hdata,"data","l");
  leg->AddEntry(hmc,"MC","l");
  leg->AddEntry(c->GetPrimitive("cfdata_py_0.0_1.0"),"0.0 #leq |#eta| < 1.0","lp");
  leg->AddEntry(c->GetPrimitive("cfdata_py_1.0_1.4"),"1.0 #leq |#eta| < 1.4","lp");
  leg->AddEntry(c->GetPrimitive("cfdata_py_1.4_1.7"),"1.4 #leq |#eta| < 1.7","lp");
  leg->AddEntry(c->GetPrimitive("cfdata_py_1.7_2.0"),"1.7 #leq |#eta| < 2.0","lp");
  leg->AddEntry(c->GetPrimitive("cfdata_py_2.0_2.5"),"2.0 #leq |#eta| < 2.5","lp");
  //leg->Draw();
  axisowner->SetTitle("2016preVFP");
  c->Draw();
  c->Update();
  ((TPaveText*)c->GetPrimitive("title"))->SetTextSize(0.075);
  c->SaveAs("fig/KPS/CFRate/ee2016a.png");

  f=new TFile("data/Run2UltraLegacy_v1/2016postVFP/SMP/CFRate_old0.root");
  c=(TCanvas*)f->Get("cf_pt");
  leg=new TLegend(0.17,0.88,0.5,0.6);
  axisowner=(TH1*)c->GetPrimitive("cfdata_py_0.0_1.0");
  axisowner->GetXaxis()->SetTitle("Electron p_{T} [GeV]");
  axisowner->GetYaxis()->SetTitle("Charg-flip rate");
  ((TH1*)c->GetPrimitive("cfdata_py_2.0_2.5"))->SetLineColor(6);
  ((TH1*)c->GetPrimitive("cfmc_py_2.0_2.5"))->SetLineColor(6);
  hdata=new TH1D;
  hmc=new TH1D;
  hmc->SetLineStyle(2);
  leg->AddEntry(hdata,"data","l");
  leg->AddEntry(hmc,"MC","l");
  leg->AddEntry(c->GetPrimitive("cfdata_py_0.0_1.0"),"0.0 #leq |#eta| < 1.0","lp");
  leg->AddEntry(c->GetPrimitive("cfdata_py_1.0_1.4"),"1.0 #leq |#eta| < 1.4","lp");
  leg->AddEntry(c->GetPrimitive("cfdata_py_1.4_1.7"),"1.4 #leq |#eta| < 1.7","lp");
  leg->AddEntry(c->GetPrimitive("cfdata_py_1.7_2.0"),"1.7 #leq |#eta| < 2.0","lp");
  leg->AddEntry(c->GetPrimitive("cfdata_py_2.0_2.5"),"2.0 #leq |#eta| < 2.5","lp");
  //leg->Draw();
  axisowner->SetTitle("2016postVFP");
  c->Draw();
  c->Update();
  ((TPaveText*)c->GetPrimitive("title"))->SetTextSize(0.075);
  c->SaveAs("fig/KPS/CFRate/ee2016b.png");
}
void plot_ref(){
  EfficiencyPlotter amc("data ^amc+tau_amc+vv+wjets+tttw");
  EfficiencyPlotter amcss("data ^amc+tau_amc+vv+wjets+tttw+ss");
  EfficiencyPlotter amcss2("data ^amc+tau_amc+vv+wjets+tttw+1.7*ss");
  EfficiencyPlotter amcS20("data ^amcS20+tau_amcS20+vv+wjets+tttw");
  EfficiencyPlotter amcS20ss("data ^amcS20+tau_amcS20+vv+wjets+tttw+ss");
  EfficiencyPlotter amcS20ss2("data ^amcS20+tau_amcS20+vv+wjets+tttw+1.7*ss");
  EfficiencyPlotter mg("data ^mg+tau_mg+vv+wjets+tttw");
  EfficiencyPlotter mgss("data ^mg+tau_mg+vv+wjets+tttw+ss");
  EfficiencyPlotter mgss2("data ^mg+tau_mg+vv+wjets+tttw+1.7*ss");
  EfficiencyPlotter mi("data ^minnlo+tau_minnlo+vv+wjets+tttw");
  EfficiencyPlotter miss("data ^minnlo+tau_minnlo+vv+wjets+tttw+ss");
  EfficiencyPlotter miss2("data ^minnlo+tau_minnlo+vv+wjets+tttw+1.7*ss");
  for(TString channel:{"mm","ee","mu","el"}){
    for(TString era:{"2016a","2016b","2017","2018"}){
      TString title=channel;
      if(era=="2016a") title+=" 2016preVFP";
      else if(era=="2016b") title+=" 2016postVFP";
      else title+=" "+era;
      TString lepton="Muon";
      if(channel[0]=='e') lepton="Electron";
      TString ll="#mu#mu";
      if(channel[0]=='e') ll="ee";
      amc.SavePlot("ref/eff/"+channel+era+"_dimass_amc","histname:"+channel+era+"/m52to150/dimass norm rebin:2 1:logy xtitle:'m("+ll+") [GeV]' title++ title:'"+title+"'");
      amc.SavePlot("ref/eff/"+channel+era+"_leta_amc","histname:"+channel+era+"/m52to150/leta norm 1:logy xmin:-2.5 xmax:2.5 noleg xtitle:'"+lepton+" #eta' title++ title:'"+title+"'");
      amc.SavePlot("ref/eff/"+channel+era+"_lpt_amc","histname:"+channel+era+"/m52to150/lpt norm 1:logy xmax:100 xtitle:'"+lepton+" p_{T} [GeV]' title++ title:'"+title+"'");
      if(channel.Contains("e")){
	amcss.SavePlot("ref/eff/"+channel+era+"_dimass_amcss","histname:"+channel+era+"/m52to150/dimass norm rebin:2 1:logy xtitle:'m("+ll+") [GeV]' title++ title:'"+title+"'");
	amcss.SavePlot("ref/eff/"+channel+era+"_leta_amcss","histname:"+channel+era+"/m52to150/leta norm 1:logy xmin:-2.5 xmax:2.5 noleg xtitle:'"+lepton+" #eta' title++ title:'"+title+"'");
	amcss.SavePlot("ref/eff/"+channel+era+"_lpt_amcss","histname:"+channel+era+"/m52to150/lpt norm 1:logy xmax:100 xtitle:'"+lepton+" p_{T} [GeV]' title++ title:'"+title+"'");
	amcss.SavePlot("ref/eff/"+channel+era+"_lsceta_amcss","histname:"+channel+era+"/m52to150/lsceta norm 1:logy xmin:-2.5 xmax:2.5 noleg xtitle:'"+lepton+" #eta_{SC}' title++ title:'"+title+"'");
	amc.SavePlot("ref/eff/"+channel+era+"_lsceta_amc","histname:"+channel+era+"/m52to150/lsceta norm 1:logy xmin:-2.5 xmax:2.5 noleg xtitle:'"+lepton+" #eta_{SC}' title++ title:'"+title+"'");
      }
      if(channel.Contains("m")){
	amcss2.SavePlot("ref/eff/"+channel+era+"_dimass_amcss2","histname:"+channel+era+"/m52to150/dimass norm rebin:2 1:logy xtitle:'m("+ll+") [GeV]' title++ title:'"+title+"'");
	amcss2.SavePlot("ref/eff/"+channel+era+"_leta_amcss2","histname:"+channel+era+"/m52to150/leta norm 1:logy xmin:-2.5 xmax:2.5 noleg xtitle:'"+lepton+" #eta' title++ title:'"+title+"'");
	amcss2.SavePlot("ref/eff/"+channel+era+"_lpt_amcss2","histname:"+channel+era+"/m52to150/lpt norm 1:logy xmax:100 xtitle:'"+lepton+" p_{T} [GeV]' title++ title:'"+title+"'");
      }

      mg.SavePlot("ref/eff/"+channel+era+"_dimass_mg","histname:"+channel+era+"/m52to150/dimass norm rebin:2 1:logy xtitle:'m("+ll+") [GeV]' title++ title:'"+title+"'");
      mg.SavePlot("ref/eff/"+channel+era+"_leta_mg","histname:"+channel+era+"/m52to150/leta norm 1:logy xmin:-2.5 xmax:2.5 noleg xtitle:'"+lepton+" #eta' title++ title:'"+title+"'");
      mg.SavePlot("ref/eff/"+channel+era+"_lpt_mg","histname:"+channel+era+"/m52to150/lpt norm 1:logy xmax:100 xtitle:'"+lepton+" p_{T} [GeV]' title++ title:'"+title+"'");
      if(channel.Contains("e")){
	mgss.SavePlot("ref/eff/"+channel+era+"_dimass_mgss","histname:"+channel+era+"/m52to150/dimass norm rebin:2 1:logy xtitle:'m("+ll+") [GeV]' title++ title:'"+title+"'");
	mgss.SavePlot("ref/eff/"+channel+era+"_leta_mgss","histname:"+channel+era+"/m52to150/leta norm 1:logy xmin:-2.5 xmax:2.5 noleg xtitle:'"+lepton+" #eta' title++ title:'"+title+"'");
	mgss.SavePlot("ref/eff/"+channel+era+"_lpt_mgss","histname:"+channel+era+"/m52to150/lpt norm 1:logy xmax:100 xtitle:'"+lepton+" p_{T} [GeV]' title++ title:'"+title+"'");
	mgss.SavePlot("ref/eff/"+channel+era+"_lsceta_mgss","histname:"+channel+era+"/m52to150/lsceta norm 1:logy xmin:-2.5 xmax:2.5 noleg xtitle:'"+lepton+" #eta_{SC}' title++ title:'"+title+"'");
	mg.SavePlot("ref/eff/"+channel+era+"_lsceta_mg","histname:"+channel+era+"/m52to150/lsceta norm 1:logy xmin:-2.5 xmax:2.5 noleg xtitle:'"+lepton+" #eta_{SC}' title++ title:'"+title+"'");
      }
      if(channel.Contains("m")){
	mgss2.SavePlot("ref/eff/"+channel+era+"_dimass_mgss2","histname:"+channel+era+"/m52to150/dimass norm rebin:2 1:logy xtitle:'m("+ll+") [GeV]' title++ title:'"+title+"'");
	mgss2.SavePlot("ref/eff/"+channel+era+"_leta_mgss2","histname:"+channel+era+"/m52to150/leta norm 1:logy xmin:-2.5 xmax:2.5 noleg xtitle:'"+lepton+" #eta' title++ title:'"+title+"'");
	mgss2.SavePlot("ref/eff/"+channel+era+"_lpt_mgss2","histname:"+channel+era+"/m52to150/lpt norm 1:logy xmax:100 xtitle:'"+lepton+" p_{T} [GeV]' title++ title:'"+title+"'");
      }

      if(!era.Contains("2016")) continue;
      mi.SavePlot("ref/eff/"+channel+era+"_dimass_mi","histname:"+channel+era+"/m52to150/dimass norm rebin:2 1:logy xtitle:'m("+ll+") [GeV]' title++ title:'"+title+"'");
      mi.SavePlot("ref/eff/"+channel+era+"_leta_mi","histname:"+channel+era+"/m52to150/leta norm 1:logy xmin:-2.5 xmax:2.5 noleg xtitle:'"+lepton+" #eta' title++ title:'"+title+"'");
      mi.SavePlot("ref/eff/"+channel+era+"_lpt_mi","histname:"+channel+era+"/m52to150/lpt norm 1:logy xmax:100 xtitle:'"+lepton+" p_{T} [GeV]' title++ title:'"+title+"'");
      if(channel.Contains("e")){
	miss.SavePlot("ref/eff/"+channel+era+"_dimass_miss","histname:"+channel+era+"/m52to150/dimass norm rebin:2 1:logy xtitle:'m("+ll+") [GeV]' title++ title:'"+title+"'");
	miss.SavePlot("ref/eff/"+channel+era+"_leta_miss","histname:"+channel+era+"/m52to150/leta norm 1:logy xmin:-2.5 xmax:2.5 noleg xtitle:'"+lepton+" #eta' title++ title:'"+title+"'");
	miss.SavePlot("ref/eff/"+channel+era+"_lpt_miss","histname:"+channel+era+"/m52to150/lpt norm 1:logy xmax:100 xtitle:'"+lepton+" p_{T} [GeV]' title++ title:'"+title+"'");
	miss.SavePlot("ref/eff/"+channel+era+"_lsceta_miss","histname:"+channel+era+"/m52to150/lsceta norm 1:logy xmin:-2.5 xmax:2.5 noleg xtitle:'"+lepton+" #eta_{SC}' title++ title:'"+title+"'");
	mi.SavePlot("ref/eff/"+channel+era+"_lsceta_mi","histname:"+channel+era+"/m52to150/lsceta norm 1:logy xmin:-2.5 xmax:2.5 noleg xtitle:'"+lepton+" #eta_{SC}' title++ title:'"+title+"'");
      }
      if(channel.Contains("m")){
	miss2.SavePlot("ref/eff/"+channel+era+"_dimass_miss2","histname:"+channel+era+"/m52to150/dimass norm rebin:2 1:logy xtitle:'m("+ll+") [GeV]' title++ title:'"+title+"'");
	miss2.SavePlot("ref/eff/"+channel+era+"_leta_miss2","histname:"+channel+era+"/m52to150/leta norm 1:logy xmin:-2.5 xmax:2.5 noleg xtitle:'"+lepton+" #eta' title++ title:'"+title+"'");
	miss2.SavePlot("ref/eff/"+channel+era+"_lpt_miss2","histname:"+channel+era+"/m52to150/lpt norm 1:logy xmax:100 xtitle:'"+lepton+" p_{T} [GeV]' title++ title:'"+title+"'");
      }
      if(!era.Contains("2016a")) continue;
      amcS20.SavePlot("ref/eff/"+channel+era+"_dimass_amcS20","histname:"+channel+era+"/m52to150/dimass norm rebin:2 1:logy xtitle:'m("+ll+") [GeV]' title++ title:'"+title+"'");
      amcS20.SavePlot("ref/eff/"+channel+era+"_leta_amcS20","histname:"+channel+era+"/m52to150/leta norm 1:logy xmin:-2.5 xmax:2.5 noleg xtitle:'"+lepton+" #eta' title++ title:'"+title+"'");
      amcS20.SavePlot("ref/eff/"+channel+era+"_lpt_amcS20","histname:"+channel+era+"/m52to150/lpt norm 1:logy xmax:100 xtitle:'"+lepton+" p_{T} [GeV]' title++ title:'"+title+"'");
      if(channel.Contains("e")){
	amcS20ss.SavePlot("ref/eff/"+channel+era+"_dimass_amcS20ss","histname:"+channel+era+"/m52to150/dimass norm rebin:2 1:logy xtitle:'m("+ll+") [GeV]' title++ title:'"+title+"'");
	amcS20ss.SavePlot("ref/eff/"+channel+era+"_leta_amcS20ss","histname:"+channel+era+"/m52to150/leta norm 1:logy xmin:-2.5 xmax:2.5 noleg xtitle:'"+lepton+" #eta' title++ title:'"+title+"'");
	amcS20ss.SavePlot("ref/eff/"+channel+era+"_lpt_amcS20ss","histname:"+channel+era+"/m52to150/lpt norm 1:logy xmax:100 xtitle:'"+lepton+" p_{T} [GeV]' title++ title:'"+title+"'");
	amcS20ss.SavePlot("ref/eff/"+channel+era+"_lsceta_amcS20ss","histname:"+channel+era+"/m52to150/lsceta norm 1:logy xmin:-2.5 xmax:2.5 noleg xtitle:'"+lepton+" #eta_{SC}' title++ title:'"+title+"'");
	amcS20.SavePlot("ref/eff/"+channel+era+"_lsceta_amcS20","histname:"+channel+era+"/m52to150/lsceta norm 1:logy xmin:-2.5 xmax:2.5 noleg xtitle:'"+lepton+" #eta_{SC}' title++ title:'"+title+"'");
      }
      if(channel.Contains("m")){
	amcS20ss2.SavePlot("ref/eff/"+channel+era+"_dimass_amcS20ss2","histname:"+channel+era+"/m52to150/dimass norm rebin:2 1:logy xtitle:'m("+ll+") [GeV]' title++ title:'"+title+"'");
	amcS20ss2.SavePlot("ref/eff/"+channel+era+"_leta_amcS20ss2","histname:"+channel+era+"/m52to150/leta norm 1:logy xmin:-2.5 xmax:2.5 noleg xtitle:'"+lepton+" #eta' title++ title:'"+title+"'");
	amcS20ss2.SavePlot("ref/eff/"+channel+era+"_lpt_amcS20ss2","histname:"+channel+era+"/m52to150/lpt norm 1:logy xmax:100 xtitle:'"+lepton+" p_{T} [GeV]' title++ title:'"+title+"'");
      }
    }
  }
}

#include"DZPlotter.cc"
TH1* GetEff(TH1* num, TH1* den, TH2* href=NULL, int iter=10){
  double numtotal,dentotal;
  if(num->InheritsFrom("TH2")){
    numtotal=((TH2*)num)->Integral(0,num->GetNbinsX()+1,0,num->GetNbinsY()+1);
    dentotal=((TH2*)den)->Integral(0,den->GetNbinsX()+1,0,den->GetNbinsY()+1);
  }else{
    numtotal=num->Integral(0,num->GetNbinsX()+1);
    dentotal=den->Integral(0,den->GetNbinsX()+1);
  }
  double mean=sqrt(numtotal/dentotal);
  cout<<"mean:"<<mean<<endl;
  TH1* heff=(TH1*)num->Clone("eff");
  heff->Divide(num,den,1,1,"B");
  for(int i=0,n=heff->GetNcells();i<n;i++){
    double val=heff->GetBinContent(i);
    double err=heff->GetBinError(i);
    if(val==0) continue;
    val/=mean;
    err/=mean;
    if(val>1.0){
      cout<<"bin "<<i<<" val: "<<val<<endl;
      val=1.0;
    }
    heff->SetBinContent(i,val);
    heff->SetBinError(i,err);
  }

  if(href){
    TH1* horigin=(TH1*)num->Clone("origin");
    horigin->Divide(num,den,1,1,"B");
    for(int it=0;it<iter;it++){
      for(int i=0,n=heff->GetNcells();i<n;i++){
	double val=horigin->GetBinContent(i);
	double err=horigin->GetBinError(i);
	if(val==0) continue;
	double sumw=0;
	double sumwx=0;
	for(int j=0;j<href->GetNbinsY()+2;j++){
	  double w=href->GetBinContent(i,j);
	  double eff=heff->GetBinContent(j);
	  if(w&&eff){
	    sumw+=w;
	    sumwx+=w*eff;
	  }
	}
	mean=sumwx/sumw;
	heff->SetBinContent(i,val/mean);
	heff->SetBinError(i,err/mean);
      }
      //cout<<"iter:"<<it<<endl;
      //for(int i=0,n=heff->GetNcells();i<n;i++){
      //cout<<heff->GetBinContent(i)<<"+"<<heff->GetBinError(i)<<" ";
      //}
      //cout<<endl;
    }
  }
  
  return heff;
}
TH2D* GetProduct2D(TH1D* hist1,TH1D* hist2){
  int nbin=hist1->GetNbinsX();
  vector<double> bins;
  for(int i=1;i<nbin+2;i++) bins.push_back(hist1->GetBinLowEdge(i));
  TH2D* rt=new TH2D("square","square",nbin,&bins[0],nbin,&bins[0]);
  for(int i=0;i<nbin+2;i++){
    for(int j=0;j<nbin+2;j++){
      double val1=hist1->GetBinContent(i);
      double err1=hist1->GetBinError(i);
      double val2=hist2->GetBinContent(j);
      double err2=hist2->GetBinError(j);
      rt->SetBinContent(i,j,val1*val2);
      if(i==j&&val1==val2&&err1&&err2)
	rt->SetBinError(i,j,val1*err2+val2*err1);
      else
	rt->SetBinError(i,j,sqrt(pow(val1*err2,2)+pow(val2*err1,2)));
    }
  }
  return rt;
}
TH2D* GetSquare2D(TH1D* hist){
  return GetProduct2D(hist,hist);
}

void Compare(TH2D* hist1,TH2D* hist2){
  for(int i=1;i<hist1->GetNbinsX()+1;i++){
    TH1D* py1=hist1->ProjectionY(Form("py1_%d",i),i,i);
    TH1D* py2=hist2->ProjectionY(Form("py2_%d",i),i,i);
    new TCanvas;
    py1->SetStats(0);
    py1->Draw("hist e");
    py2->Draw("hist e same");
    py2->SetLineColor(2);
    py2->SetMarkerStyle(0);
    TLegend *leg=new TLegend(0.9,0.9,0.7,0.8);
    leg->AddEntry(py1,hist1->GetTitle());
    leg->AddEntry(py2,hist2->GetTitle());
    leg->Draw();
    py1->SetTitle(Form("[%.1f,%.1f]",hist1->GetXaxis()->GetBinLowEdge(i),hist1->GetXaxis()->GetBinLowEdge(i+1)));
    py1->GetYaxis()->SetRangeUser(0.95,1.01);
  }
}

void SaveEff2D(TString era="ee2016a"){
  vector<TH1*> hists;
  DZPlotter aa("data mi");
  TH2D* den_data=(TH2D*)aa.GetHist(0,era+"/etaptfine_den","noproject");
  TH2D* num_data=(TH2D*)aa.GetHist(0,era+"/etaptfine_num","noproject");
  TH2D* eff_data=(TH2D*)GetEff(num_data,den_data);
  eff_data->SetNameTitle("data","eff");
  hists.push_back(eff_data);

  TH2D* den_sim=(TH2D*)aa.GetHist(1,era+"/etaptfine_den","noproject");
  TH2D* num_sim=(TH2D*)aa.GetHist(1,era+"/etaptfine_num","noproject");
  TH2D* eff_sim=(TH2D*)GetEff(num_sim,den_sim);
  eff_sim->SetNameTitle("sim","eff");
  hists.push_back(eff_sim);
  
  TH2D* sf=(TH2D*)eff_data->Clone("sf");
  sf->Divide(eff_sim);
  sf->SetNameTitle("sf","sf");
  hists.push_back(sf);

  TH4D* den_data_4Dfine=(TH4D*)aa.GetHist(0,era+"/epepfine_den","noproject");
  TH4D* num_data_4Dfine=(TH4D*)aa.GetHist(0,era+"/epepfine_num","noproject");
  TH4D* eff_data_4Dfine=(TH4D*)num_data_4Dfine->Clone("eff");
  eff_data_4Dfine->Divide(num_data_4Dfine,den_data_4Dfine,1,1,"B");
  eff_data_4Dfine->SetNameTitle("data_4Dfine","eff");
  hists.push_back(eff_data_4Dfine);
  
  TH4D* den_sim_4Dfine=(TH4D*)aa.GetHist(1,era+"/epepfine_den","noproject");
  TH4D* num_sim_4Dfine=(TH4D*)aa.GetHist(1,era+"/epepfine_num","noproject");
  TH4D* eff_sim_4Dfine=(TH4D*)num_sim_4Dfine->Clone("eff");
  eff_sim_4Dfine->Divide(num_sim_4Dfine,den_sim_4Dfine,1,1,"B");
  eff_sim_4Dfine->SetNameTitle("sim_4Dfine","eff");
  hists.push_back(eff_sim_4Dfine);
  
  TH4D* sf_4Dfine=(TH4D*)eff_data_4Dfine->Clone("sf");
  sf_4Dfine->Divide(eff_sim_4Dfine);
  sf_4Dfine->SetNameTitle("sf_4Dfine","sf");
  hists.push_back(sf_4Dfine);
  
  TH4D* den_data_4D=(TH4D*)aa.GetHist(0,era+"/epep_den","noproject");
  TH4D* num_data_4D=(TH4D*)aa.GetHist(0,era+"/epep_num","noproject");
  TH4D* eff_data_4D=(TH4D*)num_data_4D->Clone("eff");
  eff_data_4D->Divide(num_data_4D,den_data_4D,1,1,"B");
  eff_data_4D->SetNameTitle("data_4D","eff");
  hists.push_back(eff_data_4D);

  TH4D* den_sim_4D=(TH4D*)aa.GetHist(1,era+"/epep_den","noproject");
  TH4D* num_sim_4D=(TH4D*)aa.GetHist(1,era+"/epep_num","noproject");
  TH4D* eff_sim_4D=(TH4D*)num_sim_4D->Clone("eff");
  eff_sim_4D->Divide(num_sim_4D,den_sim_4D,1,1,"B");
  eff_sim_4D->SetNameTitle("sim_4D","eff");
  hists.push_back(eff_sim_4D);
  
  TH4D* sf_4D=(TH4D*)eff_data_4D->Clone("sf");
  sf_4D->Divide(eff_sim_4D);
  sf_4D->SetNameTitle("sf_4D","sf");
  hists.push_back(sf_4D);
  
  TFile f(era+"_DZ.root","recreate");
  for(auto hist:hists){
    if(hist->InheritsFrom("TH4D")) continue;
    if(hist->InheritsFrom("TH2")){
      new TCanvas;
      hist->SetOption("colz");
      hist->SetMinimum(hist->GetMinimum(0.5));
      hist->Draw();
      gPad->SetLogy();
    }
    hist->Write();
  }
}
void SaveEffDZ(TString era="ee2016a"){
  vector<TH1*> hists;
  DZPlotter aa("data mi");
  TH1D* eff_data=(TH1D*)aa.GetHist(0,era+"/lldz_eff");
  eff_data->SetNameTitle("data","eff");
  TH1D* eff_sim=(TH1D*)aa.GetHist(1,era+"/lldz_eff");
  eff_sim->SetNameTitle("sim","eff");
  TH1D* sf=(TH1D*)eff_data->Clone("sf");
  sf->Divide(eff_sim);
  sf->SetNameTitle("sf","sf");

  TFile f(era+"_DZ_DZ.root","recreate");
  eff_data->Write();
  eff_sim->Write();
  sf->Write();
  f.Close();
}  
void plot_lepton(){
  DZPlotter aa("data mi");
  aa.SavePlot("ee2016a/leta_eff","xmax:2.5 xmin:-2.5 ymin:0.94 ymax:1.01 save:ee2016a_leta_eff xtitle:'electron #eta' preliminary 1:ytitle:'DZ efficiency' 2:ytitle:'SF'");
  aa.SavePlot("ee2016b/leta_eff","xmax:2.5 xmin:-2.5 ymin:0.94 ymax:1.01 save:ee2016b_leta_eff xtitle:'electron #eta' preliminary 1:ytitle:'DZ efficiency' 2:ytitle:'SF'");
  aa.SavePlot("ee2017/leta_eff","xmax:2.5 xmin:-2.5 ymin:0.94 ymax:1.01 save:ee2017_leta_eff xtitle:'electron #eta' preliminary 1:ytitle:'DZ efficiency' 2:ytitle:'SF'");
  aa.SavePlot("ee2018/leta_eff","xmax:2.5 xmin:-2.5 ymin:0.94 ymax:1.01 save:ee2018_leta_eff xtitle:'electron #eta' preliminary 1:ytitle:'DZ efficiency' 2:ytitle:'SF'");

  aa.SavePlot("mm2016a/leta_eff","ymin:0.8 ymax:1.01 save:mm2016a_leta_eff xtitle:'muon #eta' preliminary 1:ytitle:'DZ efficiency' 2:ytitle:'SF'");
  aa.SavePlot("mm2016b/leta_eff","ymin:0.94 ymax:1.01 save:mm2016b_leta_eff xtitle:'muon #eta' preliminary 1:ytitle:'DZ efficiency' 2:ytitle:'SF'");
  aa.SavePlot("mm2017/leta_eff","ymin:0.94 ymax:1.01 save:mm2017_leta_eff xtitle:'muon #eta' preliminary 1:ytitle:'DZ efficiency' 2:ytitle:'SF'");
  aa.SavePlot("mm2018/leta_eff","ymin:0.94 ymax:1.01 save:mm2018_leta_eff xtitle:'muon #eta' preliminary 1:ytitle:'DZ efficiency' 2:ytitle:'SF'");

  aa.SavePlot("ee2016a/lpt_eff","xmax:200 xmin:10 logx ymin:0.94 ymax:1.01 save:ee2016a_lpt_eff xtitle:'electron p_{T} [GeV]' preliminary 1:ytitle:'DZ efficiency' 2:ytitle:'SF'");
  aa.SavePlot("ee2016b/lpt_eff","xmax:200 xmin:10 logx ymin:0.94 ymax:1.01 save:ee2016b_lpt_eff xtitle:'electron p_{T} [GeV]' preliminary 1:ytitle:'DZ efficiency' 2:ytitle:'SF'");
  aa.SavePlot("ee2017/lpt_eff","xmax:200 xmin:10 logx ymin:0.94 ymax:1.01 save:ee2017_lpt_eff xtitle:'electron p_{T} [GeV]' preliminary 1:ytitle:'DZ efficiency' 2:ytitle:'SF'");
  aa.SavePlot("ee2018/lpt_eff","xmax:200 xmin:10 logx ymin:0.94 ymax:1.01 save:ee2018_lpt_eff xtitle:'electron p_{T} [GeV]' preliminary 1:ytitle:'DZ efficiency' 2:ytitle:'SF'");

  aa.SavePlot("mm2016a/lpt_eff","xmax:200 xmin:10 logx ymin:0.8 ymax:1.01 save:mm2016a_lpt_eff xtitle:'muon p_{T} [GeV]' preliminary 1:ytitle:'DZ efficiency' 2:ytitle:'SF'");
  aa.SavePlot("mm2016b/lpt_eff","xmax:200 xmin:10 logx ymin:0.94 ymax:1.01 save:mm2016b_lpt_eff xtitle:'muon p_{T} [GeV]' preliminary 1:ytitle:'DZ efficiency' 2:ytitle:'SF'");
  aa.SavePlot("mm2017/lpt_eff","xmax:200 xmin:10 logx ymin:0.94 ymax:1.01 save:mm2017_lpt_eff xtitle:'muon p_{T} [GeV]' preliminary 1:ytitle:'DZ efficiency' 2:ytitle:'SF'");
  aa.SavePlot("mm2018/lpt_eff","xmax:200 xmin:10 logx ymin:0.94 ymax:1.01 save:mm2018_lpt_eff xtitle:'muon p_{T} [GeV]' preliminary 1:ytitle:'DZ efficiency' 2:ytitle:'SF'");

  aa.SavePlot("ee2016a/leta_eff_dzsf","xmax:2.5 xmin:-2.5 ymin:0.94 ymax:1.01 save:ee2016a_leta_eff_dzsf xtitle:'electron #eta' preliminary 1:ytitle:'DZ efficiency' 2:ytitle:'SF'");
  aa.SavePlot("ee2016b/leta_eff_dzsf","xmax:2.5 xmin:-2.5 ymin:0.94 ymax:1.01 save:ee2016b_leta_eff_dzsf xtitle:'electron #eta' preliminary 1:ytitle:'DZ efficiency' 2:ytitle:'SF'");
  aa.SavePlot("ee2017/leta_eff_dzsf","xmax:2.5 xmin:-2.5 ymin:0.94 ymax:1.01 save:ee2017_leta_eff_dzsf xtitle:'electron #eta' preliminary 1:ytitle:'DZ efficiency' 2:ytitle:'SF'");
  aa.SavePlot("ee2018/leta_eff_dzsf","xmax:2.5 xmin:-2.5 ymin:0.94 ymax:1.01 save:ee2018_leta_eff_dzsf xtitle:'electron #eta' preliminary 1:ytitle:'DZ efficiency' 2:ytitle:'SF'");

  aa.SavePlot("mm2016a/leta_eff_dzsf","ymin:0.8 ymax:1.01 save:mm2016a_leta_eff_dzsf xtitle:'muon #eta' preliminary 1:ytitle:'DZ efficiency' 2:ytitle:'SF'");
  aa.SavePlot("mm2016b/leta_eff_dzsf","ymin:0.94 ymax:1.01 save:mm2016b_leta_eff_dzsf xtitle:'muon #eta' preliminary 1:ytitle:'DZ efficiency' 2:ytitle:'SF'");
  aa.SavePlot("mm2017/leta_eff_dzsf","ymin:0.94 ymax:1.01 save:mm2017_leta_eff_dzsf xtitle:'muon #eta' preliminary 1:ytitle:'DZ efficiency' 2:ytitle:'SF'");
  aa.SavePlot("mm2018/leta_eff_dzsf","ymin:0.94 ymax:1.01 save:mm2018_leta_eff_dzsf xtitle:'muon #eta' preliminary 1:ytitle:'DZ efficiency' 2:ytitle:'SF'");

  aa.SavePlot("ee2016a/lpt_eff_dzsf","xmax:200 xmin:10 logx ymin:0.94 ymax:1.01 save:ee2016a_lpt_eff_dzsf xtitle:'electron p_{T} [GeV]' preliminary 1:ytitle:'DZ efficiency' 2:ytitle:'SF'");
  aa.SavePlot("ee2016b/lpt_eff_dzsf","xmax:200 xmin:10 logx ymin:0.94 ymax:1.01 save:ee2016b_lpt_eff_dzsf xtitle:'electron p_{T} [GeV]' preliminary 1:ytitle:'DZ efficiency' 2:ytitle:'SF'");
  aa.SavePlot("ee2017/lpt_eff_dzsf","xmax:200 xmin:10 logx ymin:0.94 ymax:1.01 save:ee2017_lpt_eff_dzsf xtitle:'electron p_{T} [GeV]' preliminary 1:ytitle:'DZ efficiency' 2:ytitle:'SF'");
  aa.SavePlot("ee2018/lpt_eff_dzsf","xmax:200 xmin:10 logx ymin:0.94 ymax:1.01 save:ee2018_lpt_eff_dzsf xtitle:'electron p_{T} [GeV]' preliminary 1:ytitle:'DZ efficiency' 2:ytitle:'SF'");

  aa.SavePlot("mm2016a/lpt_eff_dzsf","xmax:200 xmin:10 logx ymin:0.8 ymax:1.01 save:mm2016a_lpt_eff_dzsf xtitle:'muon p_{T} [GeV]' preliminary 1:ytitle:'DZ efficiency' 2:ytitle:'SF'");
  aa.SavePlot("mm2016b/lpt_eff_dzsf","xmax:200 xmin:10 logx ymin:0.94 ymax:1.01 save:mm2016b_lpt_eff_dzsf xtitle:'muon p_{T} [GeV]' preliminary 1:ytitle:'DZ efficiency' 2:ytitle:'SF'");
  aa.SavePlot("mm2017/lpt_eff_dzsf","xmax:200 xmin:10 logx ymin:0.94 ymax:1.01 save:mm2017_lpt_eff_dzsf xtitle:'muon p_{T} [GeV]' preliminary 1:ytitle:'DZ efficiency' 2:ytitle:'SF'");
  aa.SavePlot("mm2018/lpt_eff_dzsf","xmax:200 xmin:10 logx ymin:0.94 ymax:1.01 save:mm2018_lpt_eff_dzsf xtitle:'muon p_{T} [GeV]' preliminary 1:ytitle:'DZ efficiency' 2:ytitle:'SF'");
}
void plot_dirap(){
  DZPlotter aa("data mi");
  aa.SavePlot("ee2016a/dirap_eff","xmax:2.4 xmin:-2.4 ymin:0.91 ymax:1.09 save:ee2016a_dirap_eff xtitle:'y(ee)' preliminary 1:ytitle:'DZ efficiency' 2:ytitle:'SF'");
  aa.SavePlot("ee2016b/dirap_eff","xmax:2.4 xmin:-2.4 ymin:0.91 ymax:1.09 save:ee2016b_dirap_eff xtitle:'y(ee)' preliminary 1:ytitle:'DZ efficiency' 2:ytitle:'SF'");
  aa.SavePlot("ee2017/dirap_eff","xmax:2.4 xmin:-2.4 ymin:0.91 ymax:1.09 save:ee2017_dirap_eff xtitle:'y(ee)' preliminary 1:ytitle:'DZ efficiency' 2:ytitle:'SF'");
  aa.SavePlot("ee2018/dirap_eff","xmax:2.4 xmin:-2.4 ymin:0.91 ymax:1.09 save:ee2018_dirap_eff xtitle:'y(ee)' preliminary 1:ytitle:'DZ efficiency' 2:ytitle:'SF'");

  aa.SavePlot("mm2016a/dirap_eff","xmax:2.4 xmin:-2.4 ymin:0.8 ymax:1.09 save:mm2016a_dirap_eff xtitle:'y(#mu#mu)' preliminary 1:ytitle:'DZ efficiency' 2:ytitle:'SF'");
  aa.SavePlot("mm2016b/dirap_eff","xmax:2.4 xmin:-2.4 ymin:0.91 ymax:1.09 save:mm2016b_dirap_eff xtitle:'y(#mu#mu)' preliminary 1:ytitle:'DZ efficiency' 2:ytitle:'SF'");
  aa.SavePlot("mm2017/dirap_eff","xmax:2.4 xmin:-2.4 ymin:0.91 ymax:1.09 save:mm2017_dirap_eff xtitle:'y(#mu#mu)' preliminary 1:ytitle:'DZ efficiency' 2:ytitle:'SF'");
  aa.SavePlot("mm2018/dirap_eff","xmax:2.4 xmin:-2.4 ymin:0.91 ymax:1.09 save:mm2018_dirap_eff xtitle:'y(#mu#mu)' preliminary 1:ytitle:'DZ efficiency' 2:ytitle:'SF'");

  aa.SavePlot("ee2016a/dirap_eff_dzsf","xmax:2.4 xmin:-2.4 ymin:0.91 ymax:1.09 save:ee2016a_dirap_eff_dzsf xtitle:'y(ee)' preliminary 1:ytitle:'DZ efficiency' 2:ytitle:'SF'");
  aa.SavePlot("ee2016b/dirap_eff_dzsf","xmax:2.4 xmin:-2.4 ymin:0.91 ymax:1.09 save:ee2016b_dirap_eff_dzsf xtitle:'y(ee)' preliminary 1:ytitle:'DZ efficiency' 2:ytitle:'SF'");
  aa.SavePlot("ee2017/dirap_eff_dzsf","xmax:2.4 xmin:-2.4 ymin:0.91 ymax:1.09 save:ee2017_dirap_eff_dzsf xtitle:'y(ee)' preliminary 1:ytitle:'DZ efficiency' 2:ytitle:'SF'");
  aa.SavePlot("ee2018/dirap_eff_dzsf","xmax:2.4 xmin:-2.4 ymin:0.91 ymax:1.09 save:ee2018_dirap_eff_dzsf xtitle:'y(ee)' preliminary 1:ytitle:'DZ efficiency' 2:ytitle:'SF'");

  aa.SavePlot("mm2016a/dirap_eff_dzsf","xmax:2.4 xmin:-2.4 ymin:0.8 ymax:1.09 save:mm2016a_dirap_eff_dzsf xtitle:'y(#mu#mu)' preliminary 1:ytitle:'DZ efficiency' 2:ytitle:'SF'");
  aa.SavePlot("mm2016b/dirap_eff_dzsf","xmax:2.4 xmin:-2.4 ymin:0.91 ymax:1.09 save:mm2016b_dirap_eff_dzsf xtitle:'y(#mu#mu)' preliminary 1:ytitle:'DZ efficiency' 2:ytitle:'SF'");
  aa.SavePlot("mm2017/dirap_eff_dzsf","xmax:2.4 xmin:-2.4 ymin:0.91 ymax:1.09 save:mm2017_dirap_eff_dzsf xtitle:'y(#mu#mu)' preliminary 1:ytitle:'DZ efficiency' 2:ytitle:'SF'");
  aa.SavePlot("mm2018/dirap_eff_dzsf","xmax:2.4 xmin:-2.4 ymin:0.91 ymax:1.09 save:mm2018_dirap_eff_dzsf xtitle:'y(#mu#mu)' preliminary 1:ytitle:'DZ efficiency' 2:ytitle:'SF'");
}
  
  
void plot2(TString era="ee2016a",int iter=0){
  DZPlotter aa("data mi");
  TH1D* den=(TH1D*)aa.GetHist(0,era+"/leta_den","noproject");
  TH1D* num=(TH1D*)aa.GetHist(0,era+"/leta_num","noproject");
  TH2D* ref2d=(TH2D*)aa.GetHist(0,era+"/ptpt_den","noproject");
  TH1D* eff=(TH1D*)GetEff(num,den,ref2d,iter);
  TH2D* eff2d=GetSquare2D(eff);
  TH2D* eff2d_raw=(TH2D*)aa.GetHist(0,era+"/etaeta_eff","noproject");
  new TCanvas;
  eff2d->Draw("colz");
  eff2d->SetMaximum(1.01);
  eff2d->SetMinimum(0.9);
  new TCanvas;
  eff2d_raw->Draw("colz");
  eff2d_raw->SetMaximum(1.01);
  eff2d_raw->SetMinimum(0.9);

  TH1D* bias=new TH1D("bias","bias",100,-10,10);
  for(int i=1;i<eff2d->GetNbinsX()+1;i++){
    for(int j=1;j<eff2d->GetNbinsY()+1;j++){
      double val=eff2d->GetBinContent(i,j)-eff2d_raw->GetBinContent(i,j);
      double err=eff2d_raw->GetBinError(i,j);
      if(err) bias->Fill(val/err);
    }
  }
  new TCanvas;
  bias->Draw();
  bias->Fit("gaus");
  Compare(eff2d,eff2d_raw);
}  

void plot3(TString era="ee2016a",int iter=0){
  DZPlotter aa("data mi");

  TH2D* eff2d_raw_data=(TH2D*)aa.GetHist(0,era+"/ptpt_eff","noproject");
  TH2D* ref2d_data=(TH2D*)aa.GetHist(0,era+"/ptpt_den","noproject");
  TH1D* den_data=(TH1D*)aa.GetHist(0,era+"/lpt_den","noproject");
  TH1D* num_data=(TH1D*)aa.GetHist(0,era+"/lpt_num","noproject");
  TH1D* eff_data=(TH1D*)GetEff(num_data,den_data,ref2d_data,iter);
  TH2D* eff2d_data=GetSquare2D(eff_data);
  new TCanvas;
  eff2d_data->Draw("colz");
  eff2d_data->SetMaximum(1.01);
  eff2d_data->SetMinimum(0.96);
  new TCanvas;
  eff2d_raw_data->Draw("colz");
  eff2d_raw_data->SetMaximum(1.01);
  eff2d_raw_data->SetMinimum(0.96);

  TH1D* bias_data=new TH1D("bias_data","bias_data",100,-10,10);
  TH2D* bias2d_data=(TH2D*)eff2d_data->Clone("bias2d_data");
  bias2d_data->Reset();
  for(int i=1;i<eff2d_data->GetNbinsX()+1;i++){
    for(int j=1;j<eff2d_data->GetNbinsY()+1;j++){
      double val=eff2d_data->GetBinContent(i,j)-eff2d_raw_data->GetBinContent(i,j);
      double err=eff2d_raw_data->GetBinError(i,j);
      if(err){
	bias_data->Fill(val/err);
	bias2d_data->SetBinContent(i,j,val/err);
      }
    }
  }
  new TCanvas;
  bias_data->Draw();
  bias_data->Fit("gaus");
  new TCanvas;
  bias2d_data->Draw("colz");
  bias2d_data->SetMinimum(-5);
  bias2d_data->SetMaximum(5);

  TH2D* eff2d_raw_sim=(TH2D*)aa.GetHist(1,era+"/ptpt_eff","noproject");
  TH2D* ref2d_sim=(TH2D*)aa.GetHist(1,era+"/ptpt_den","noproject");
  TH1D* den_sim=(TH1D*)aa.GetHist(1,era+"/lpt_den","noproject");
  TH1D* num_sim=(TH1D*)aa.GetHist(1,era+"/lpt_num","noproject");
  TH1D* eff_sim=(TH1D*)GetEff(num_sim,den_sim,ref2d_sim,iter);
  TH2D* eff2d_sim=GetSquare2D(eff_sim);
  new TCanvas;
  eff2d_sim->Draw("colz");
  eff2d_sim->SetMaximum(1.01);
  eff2d_sim->SetMinimum(0.96);
  new TCanvas;
  eff2d_raw_sim->Draw("colz");
  eff2d_raw_sim->SetMaximum(1.01);
  eff2d_raw_sim->SetMinimum(0.96);

  TH1D* bias_sim=new TH1D("bias_sim","bias_sim",100,-10,10);
  TH2D* bias2d_sim=(TH2D*)eff2d_sim->Clone("bias2d_sim");
  bias2d_sim->Reset();
  for(int i=1;i<eff2d_sim->GetNbinsX()+1;i++){
    for(int j=1;j<eff2d_sim->GetNbinsY()+1;j++){
      double val=eff2d_sim->GetBinContent(i,j)-eff2d_raw_sim->GetBinContent(i,j);
      double err=eff2d_raw_sim->GetBinError(i,j);
      if(err){
	bias_sim->Fill(val/err);
	bias2d_sim->SetBinContent(i,j,val/err);
      }
    }
  }
  new TCanvas;
  bias_sim->Draw();
  bias_sim->Fit("gaus");
  new TCanvas;
  bias2d_sim->Draw("colz");
  bias2d_sim->SetMinimum(-5);
  bias2d_sim->SetMaximum(5);

  //////////////////
  TH2D* sf2d=(TH2D*)eff2d_data->Clone("sf2d");
  sf2d->Divide(eff2d_sim);
  TH2D* sf2d_raw=(TH2D*)eff2d_raw_data->Clone("sf2d_raw");
  sf2d_raw->Divide(eff2d_raw_sim);
  new TCanvas;
  sf2d->Draw("colz");
  sf2d->SetMaximum(1.01);
  sf2d->SetMinimum(0.96);
  new TCanvas;
  sf2d_raw->Draw("colz");
  sf2d_raw->SetMaximum(1.01);
  sf2d_raw->SetMinimum(0.96);

  TH1D* bias_sf=new TH1D("bias_sf","bias_sf",100,-10,10);
  TH2D* bias2d_sf=(TH2D*)sf2d->Clone("bias2d_sf");
  bias2d_sf->Reset();
  for(int i=1;i<sf2d->GetNbinsX()+1;i++){
    for(int j=1;j<sf2d->GetNbinsY()+1;j++){
      double val=sf2d->GetBinContent(i,j)-sf2d_raw->GetBinContent(i,j);
      //double err=sqrt(pow(sf2d->GetBinError(i,j),2)+pow(sf2d_raw->GetBinError(i,j),2));
      double err=sf2d_raw->GetBinError(i,j);
      if(err&&sf2d->GetBinContent(i,j)){
	bias_sf->Fill(val/err);
	bias2d_sf->SetBinContent(i,j,val/err);
      }
    }
  }
  new TCanvas;
  bias_sf->Draw();
  bias_sf->Fit("gaus");
  new TCanvas;
  bias2d_sf->Draw("colz");
  bias2d_sf->SetMinimum(-5);
  bias2d_sf->SetMaximum(5);
  
  sf2d->SetTitle("1D#times1D");
  sf2d_raw->SetTitle("2D");
  Compare(sf2d,sf2d_raw);
}  

void plot4(TString era="ee2016a",int iter=0){
  DZPlotter aa("data mi");

  TH2D* eff2d_raw_data=(TH2D*)aa.GetHist(0,era+"/etaeta_eff","noproject");
  TH2D* ref2d_data=(TH2D*)aa.GetHist(0,era+"/etaeta_den","noproject");
  TH1D* den_data=(TH1D*)aa.GetHist(0,era+"/leta_den","noproject");
  TH1D* num_data=(TH1D*)aa.GetHist(0,era+"/leta_num","noproject");
  TH1D* eff_data=(TH1D*)GetEff(num_data,den_data,ref2d_data,iter);
  TH2D* eff2d_data=GetSquare2D(eff_data);
  new TCanvas;
  eff2d_data->Draw("colz");
  eff2d_data->SetMaximum(1.01);
  eff2d_data->SetMinimum(0.96);
  new TCanvas;
  eff2d_raw_data->Draw("colz");
  eff2d_raw_data->SetMaximum(1.01);
  eff2d_raw_data->SetMinimum(0.96);

  TH1D* bias_data=new TH1D("bias_data","bias_data",100,-10,10);
  TH2D* bias2d_data=(TH2D*)eff2d_data->Clone("bias2d_data");
  bias2d_data->Reset();
  for(int i=1;i<eff2d_data->GetNbinsX()+1;i++){
    for(int j=1;j<eff2d_data->GetNbinsY()+1;j++){
      double val=eff2d_data->GetBinContent(i,j)-eff2d_raw_data->GetBinContent(i,j);
      double err=eff2d_raw_data->GetBinError(i,j);
      if(err){
	bias_data->Fill(val/err);
	bias2d_data->SetBinContent(i,j,val/err);
      }
    }
  }
  new TCanvas;
  bias_data->Draw();
  bias_data->Fit("gaus");
  new TCanvas;
  bias2d_data->Draw("colz");
  bias2d_data->SetMinimum(-5);
  bias2d_data->SetMaximum(5);

  TH2D* eff2d_raw_sim=(TH2D*)aa.GetHist(1,era+"/etaeta_eff","noproject");
  TH2D* ref2d_sim=(TH2D*)aa.GetHist(1,era+"/etaeta_den","noproject");
  TH1D* den_sim=(TH1D*)aa.GetHist(1,era+"/leta_den","noproject");
  TH1D* num_sim=(TH1D*)aa.GetHist(1,era+"/leta_num","noproject");
  TH1D* eff_sim=(TH1D*)GetEff(num_sim,den_sim,ref2d_sim,iter);
  TH2D* eff2d_sim=GetSquare2D(eff_sim);
  new TCanvas;
  eff2d_sim->Draw("colz");
  eff2d_sim->SetMaximum(1.01);
  eff2d_sim->SetMinimum(0.96);
  new TCanvas;
  eff2d_raw_sim->Draw("colz");
  eff2d_raw_sim->SetMaximum(1.01);
  eff2d_raw_sim->SetMinimum(0.96);

  TH1D* bias_sim=new TH1D("bias_sim","bias_sim",100,-10,10);
  TH2D* bias2d_sim=(TH2D*)eff2d_sim->Clone("bias2d_sim");
  bias2d_sim->Reset();
  for(int i=1;i<eff2d_sim->GetNbinsX()+1;i++){
    for(int j=1;j<eff2d_sim->GetNbinsY()+1;j++){
      double val=eff2d_sim->GetBinContent(i,j)-eff2d_raw_sim->GetBinContent(i,j);
      double err=eff2d_raw_sim->GetBinError(i,j);
      if(err){
	bias_sim->Fill(val/err);
	bias2d_sim->SetBinContent(i,j,val/err);
      }
    }
  }
  new TCanvas;
  bias_sim->Draw();
  bias_sim->Fit("gaus");
  new TCanvas;
  bias2d_sim->Draw("colz");
  bias2d_sim->SetMinimum(-5);
  bias2d_sim->SetMaximum(5);

  //////////////////
  TH2D* sf2d=(TH2D*)eff2d_data->Clone("sf2d");
  sf2d->Divide(eff2d_sim);
  TH2D* sf2d_raw=(TH2D*)eff2d_raw_data->Clone("sf2d_raw");
  sf2d_raw->Divide(eff2d_raw_sim);
  new TCanvas;
  sf2d->Draw("colz");
  sf2d->SetMaximum(1.01);
  sf2d->SetMinimum(0.96);
  new TCanvas;
  sf2d_raw->Draw("colz");
  sf2d_raw->SetMaximum(1.01);
  sf2d_raw->SetMinimum(0.96);

  TH1D* bias_sf=new TH1D("bias_sf","bias_sf",100,-10,10);
  TH2D* bias2d_sf=(TH2D*)sf2d->Clone("bias2d_sf");
  bias2d_sf->Reset();
  for(int i=1;i<sf2d->GetNbinsX()+1;i++){
    for(int j=1;j<sf2d->GetNbinsY()+1;j++){
      double val=sf2d->GetBinContent(i,j)-sf2d_raw->GetBinContent(i,j);
      //double err=sqrt(pow(sf2d->GetBinError(i,j),2)+pow(sf2d_raw->GetBinError(i,j),2));
      double err=sf2d_raw->GetBinError(i,j);
      if(err&&sf2d->GetBinContent(i,j)){
	bias_sf->Fill(val/err);
	bias2d_sf->SetBinContent(i,j,val/err);
      }
    }
  }
  new TCanvas;
  bias_sf->Draw();
  bias_sf->SetTitle("sf pull, "+era);
  bias_sf->Fit("gaus");
  new TCanvas;
  bias2d_sf->Draw("colz");
  bias2d_sf->SetMinimum(-5);
  bias2d_sf->SetMaximum(5);
  
  eff2d_data->SetTitle("1D#times1D");
  eff2d_raw_data->SetTitle("2D");
  Compare(eff2d_data,eff2d_raw_data);

  eff2d_sim->SetTitle("1D#times1D");
  eff2d_raw_sim->SetTitle("2D");
  //Compare(eff2d_sim,eff2d_raw_sim);

  sf2d->SetTitle("1D#times1D");
  sf2d_raw->SetTitle("2D");
  //Compare(sf2d,sf2d_raw);
}  

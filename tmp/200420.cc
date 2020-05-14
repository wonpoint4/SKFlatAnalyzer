void validation(){
  EfficiencyPlotter a;
  vector<TString> plots={"l0pt","l1pt","leta","lsceta","dirap","dimass"};
  map<TString,vector<TString>> channels={
    {"ee2016",{"MediumID_Q"}},
    {"ee2017",{"MediumID_Q"}},
    {"ee2018",{"MediumID_Q"}},
    {"el2016",{"MediumID_Q","TightID_Selective_Q"}},
    {"el2017",{"MediumID_Q","TightID_Selective_Q"}},
    {"el2018",{"MediumID_Q","TightID_Selective_Q"}},
    {"mm2016",{"MediumID_trkIsoLoose_Q"}},
    {"mm2017",{"MediumID_trkIsoLoose_Q"}},
    {"mm2018",{"MediumID_trkIsoLoose_Q"}},
    {"mu2016",{"MediumID_trkIsoLoose_Q"}},
    {"mu2017",{"MediumID_trkIsoLoose_Q"}},
    {"mu2018",{"MediumID_trkIsoLoose_Q"}},
  };
        
  for(auto [channel,ids]:channels){
    for(auto id:ids){
      for(auto plot:plots){
	if(channel.Contains(TRegexp("^m"))&&plot=="lsceta") continue;
	TString option="norm sysname:efficiencySF";
	if(plot.Contains("pt")) option+=" xmax:70";
	a.SavePlot(channel+"/m80to100/"+plot+"_"+id,option);
      }
    }
  }
}
void normalization(){
  EfficiencyPlotter a;
  map<TString,vector<TString>> channels={
    {"ee2016",{"MediumID_Q"}},
    {"ee2017",{"MediumID_Q"}},
    {"ee2018",{"MediumID_Q"}},
    {"el2016",{"MediumID_Q","TightID_Selective_Q"}},
    {"el2017",{"MediumID_Q","TightID_Selective_Q"}},
    {"el2018",{"MediumID_Q","TightID_Selective_Q"}},
    {"mm2016",{"MediumID_trkIsoLoose_Q"}},
    {"mm2017",{"MediumID_trkIsoLoose_Q"}},
    {"mm2018",{"MediumID_trkIsoLoose_Q"}},
    {"mu2016",{"MediumID_trkIsoLoose_Q"}},
    {"mu2017",{"MediumID_trkIsoLoose_Q"}},
    {"mu2018",{"MediumID_trkIsoLoose_Q"}},
  };
  DEBUG=0;
  for(auto [channel,ids]:channels){
    for(auto id:ids){
      TCanvas* c=a.DrawPlot(channel+"/m80to100/dimass_"+id,"type:1");
      TH1D* hdata=(TH1D*)c->FindObject("data");
      TH1D* hsim=(TH1D*)c->FindObject("sim");
      if(strcmp(hsim->ClassName(),"TH1D")!=0){
	hsim->SetName("sim_stack");
	hsim=(TH1D*)c->FindObject("sim");
      }
      double ndata=hdata->Integral();
      double nsim=hsim->Integral();
      
      cout<<channel<<" "<<id<<" "<<ndata<<"\t"<<nsim<<"\t"<<ndata/nsim<<endl;
    }
  }
}
    

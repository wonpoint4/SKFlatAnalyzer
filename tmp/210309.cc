{
  EfficiencyPlotter a;
  TString prefixes[]={"mm2016a","mm2016b","mm2017","mm2018",
		  "mu2016a","mu2016b","mu2017","mu2018",
		  "ee2016a","ee2016b","ee2017","ee2018",
		  "el2016a","el2016b","el2017","el2018",};
  for(auto prefix:prefixes){
    a.SavePlot("plot/"+prefix+"_lpt","histname:"+prefix+"/m80to100/lpt xmax:80 norm xtitle:'lepton p_{T}'");
    a.SavePlot("plot/"+prefix+"_l0pt","histname:"+prefix+"/m80to100/l0pt xmax:80 norm xtitle:'leading lepton p_{T}'");
    a.SavePlot("plot/"+prefix+"_l1pt","histname:"+prefix+"/m80to100/l1pt xmax:80 norm xtitle:'sub-leading lepton p_{T}'");
    a.SavePlot("plot/"+prefix+"_leta","histname:"+prefix+"/m80to100/leta norm xmin:-2.5 xmax:2.5 noleg xtitle:'lepton #eta'");
    if(prefix.Contains("el")||prefix.Contains("ee")) a.SavePlot("plot/"+prefix+"_lsceta","histname:"+prefix+"/m80to100/lsceta norm xmin:-2.5 xmax:2.5 noleg xtitle:'lepton #eta_{SC}'");
    a.SavePlot("plot/"+prefix+"_dipt","histname:"+prefix+"/m80to100/dipt norm xmax:200 rebin:2 xtitle:p_{T}(ll)");
  }
  EfficiencyPlotter b("data ^amc+tau_amc+vv+wjets+tttw+2*ss");
  for(auto prefix:prefixes){
    if(prefix.Contains("el")||prefix.Contains("ee")) continue;
    b.SavePlot("plot/"+prefix+"_dimass","histname:"+prefix+"/m60to120/dimass norm xtitle:m(ll)");
  }

  EfficiencyPlotter c("data ^amc+tau_amc+vv+wjets+tttw+ss");
  for(auto prefix:prefixes){
    if(prefix.Contains("mu")||prefix.Contains("mm")) continue;
    c.SavePlot("plot/"+prefix+"_dimass","histname:"+prefix+"/m60to120/dimass norm xtitle:m(ll)");
  }  
}

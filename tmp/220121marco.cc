{
  EfficiencyPlotter aa("data ^mi+tau_mi+vv+wjets+tttw+1.7*ss_mi");
  aa.SavePlot("marco/mu2016a_m80to100_dirap","histname:mu2016a/m80to100/dirap norm rebin:2 2:ymin:0.901 2:ymax:1.099 BMleg preliminary xtitle:'y_{#mu#mu}'");
  aa.SavePlot("marco/mu2016a_m80to100_l0eta","histname:mu2016a/m80to100/l0eta norm rebin:2 2:ymin:0.901 2:ymax:1.099 BMleg preliminary xtitle:'Leading #mu p_{T} [GeV]'");
  aa.SavePlot("marco/mu2016a_m80to100_l1eta","histname:mu2016a/m80to100/l1eta norm rebin:2 2:ymin:0.901 2:ymax:1.099 BMleg preliminary xtitle:'Sub-leading #mu p_{T} [GeV]'");
  aa.SavePlot("marco/mu2016b_m80to100_dirap","histname:mu2016b/m80to100/dirap norm rebin:2 2:ymin:0.901 2:ymax:1.099 BMleg preliminary xtitle:'y_{#mu#mu}'");
  aa.SavePlot("marco/mu2016b_m80to100_l0eta","histname:mu2016b/m80to100/l0eta norm rebin:2 2:ymin:0.901 2:ymax:1.099 BMleg preliminary xtitle:'Leading #mu p_{T} [GeV]'");
  aa.SavePlot("marco/mu2016b_m80to100_l1eta","histname:mu2016b/m80to100/l1eta norm rebin:2 2:ymin:0.901 2:ymax:1.099 BMleg preliminary xtitle:'Sub-leading #mu p_{T} [GeV]'");
}

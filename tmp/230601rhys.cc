#include"EfficiencyPlotter.cc"
void plot(){
  EfficiencyPlotter aa;
  aa.plotdir="fig/230601rhys";
  aa.SavePlot("ee2017/m80to100/leta","norm xmin:-2.5 xmax:2.5 xtitle:'electron #eta' save:ee2017_leta_central");
  aa.SavePlot("ee2017/m80to100/leta","norm xmin:-2.5 xmax:2.5 suffix:_newprefireweight3:sim xtitle:'electron #eta' save:ee2017_leta_mode3");
  aa.SavePlot("ee2017/m80to100/leta","norm xmin:-2.5 xmax:2.5 suffix:_newprefireweight4:sim xtitle:'electron #eta' save:ee2017_leta_mode4");
  aa.SavePlot("ee2017/m80to100/leta","norm xmin:-2.5 xmax:2.5 suffix:_newprefireweight5:sim xtitle:'electron #eta' save:ee2017_leta_mode5");
  aa.SavePlot("ee2017/m80to100/leta","norm xmin:-2.5 xmax:2.5 suffix:_newprefireweight6:sim xtitle:'electron #eta' save:ee2017_leta_mode6");
  aa.SavePlot("ee2017/m80to100/lpt_eta2","norm xtitle:'electron p_{T} [GeV]' save:ee2017_lpt_eta2_central 1:text:0.16,0.8,'#eta>2.0'");
  aa.SavePlot("ee2017/m80to100/lpt_eta2","norm suffix:_newprefireweight3:sim xtitle:'electron p_{T} [GeV]' save:ee2017_lpt_eta2_mode3 1:text:0.16,0.8,'#eta>2.0'");
  aa.SavePlot("ee2017/m80to100/lpt_eta2","norm suffix:_newprefireweight4:sim xtitle:'electron p_{T} [GeV]' save:ee2017_lpt_eta2_mode4 1:text:0.16,0.8,'#eta>2.0'");
  aa.SavePlot("ee2017/m80to100/lpt_eta2","norm suffix:_newprefireweight5:sim xtitle:'electron p_{T} [GeV]' save:ee2017_lpt_eta2_mode5 1:text:0.16,0.8,'#eta>2.0'");
  aa.SavePlot("ee2017/m80to100/lpt_eta2","norm suffix:_newprefireweight6:sim xtitle:'electron p_{T} [GeV]' save:ee2017_lpt_eta2_mode6 1:text:0.16,0.8,'#eta>2.0'");
}

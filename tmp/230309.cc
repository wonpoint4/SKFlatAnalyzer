#include"EfficiencyPlotter.cc"

void plot(){
  EfficiencyPlotter aa;
  aa.plotdir="fig/230309";

  aa.SavePlot("el2018/m80to100/l0pt","norm xmax:100 save:el2018_l0pt xtitle:'Leading electron p_{T} [GeV]'");
  aa.SavePlot("el2018/m80to100/l0pt_nocor","norm xmax:100 save:el2018_l0pt_nocor xtitle:'Leading electron p_{T} [GeV]'");
  aa.SavePlot("el2018/m80to100/l0pt","norm xmax:100 save:el2018_l0pt_newSF suffix:_newtriggerSF:sim xtitle:'Leading electron p_{T} [GeV]'");

  aa.SavePlot("el201832/m80to100/l0pt","norm xmax:100 save:el201832_l0pt xtitle:'Leading electron p_{T} [GeV]'");
  aa.SavePlot("el201832/m80to100/l0pt_nocor","norm xmax:100 save:el201832_l0pt_nocor xtitle:'Leading electron p_{T} [GeV]'");
  aa.SavePlot("el201832/m80to100/l0pt","norm xmax:100 save:el201832_l0pt_newSF suffix:_newtriggerSF:sim xtitle:'Leading electron p_{T} [GeV]'");

  aa.SavePlot("el201828/m80to100/l0pt","norm xmax:100 save:el201828_l0pt xtitle:'Leading electron p_{T} [GeV]'");
  aa.SavePlot("el201828/m80to100/l0pt_nocor","norm xmax:100 save:el201828_l0pt_nocor xtitle:'Leading electron p_{T} [GeV]'");
  aa.SavePlot("el201828/m80to100/l0pt","norm xmax:100 save:el201828_l0pt_newSF suffix:_newtriggerSF:sim xtitle:'Leading electron p_{T} [GeV]'");
  
}

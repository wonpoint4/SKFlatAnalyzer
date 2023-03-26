//Electron v12 and v13 transition. with electron energy correction in eff measurement and try other sf methods

#include"EfficiencyPlotter.cc"

void plot(){
  EfficiencyPlotter aa;
  aa.plotdir="fig/230312/";
  aa.SavePlot("el2018/m80to100/l0pt","norm xmax:100 save:el2018_l0pt.png");
  aa.SavePlot("el2018/m80to100/l0pt","norm xmax:100 suffix:_newtriggerSF:dy save:el2018_l0pt_aepcor.png");
  aa.SavePlot("el2018/m80to100/l0pt","norm xmax:100 suffix:_newtriggerSF_t1:dy save:el2018_l0pt_aepcor_altSF.png");
  aa.SavePlot("el2018/m80to100/l0pt","norm xmax:100 suffix:_newtriggerSF_t2:dy save:el2018_l0pt_aepcor_Ele32N.png");
  aa.SavePlot("el2018/m80to100/l0pt","norm xmax:100 suffix:_newtriggerSF_s3:dy save:el2018_l0pt_aepcor_tagpt45.png");

  aa.SavePlot("el201828/m80to100/l0pt","norm xmax:100 save:el201828_l0pt.png");
  aa.SavePlot("el201828/m80to100/l0pt","norm xmax:100 suffix:_newtriggerSF:dy save:el201828_l0pt_aepcor.png");
  aa.SavePlot("el201828/m80to100/l0pt","norm xmax:100 suffix:_newtriggerSF_s3:dy save:el201828_l0pt_aepcor_tagpt45.png");

  aa.SavePlot("el201832/m80to100/l0pt","norm xmax:100 save:el201832_l0pt.png");
  aa.SavePlot("el201832/m80to100/l0pt","norm xmax:100 suffix:_newtriggerSF:dy save:el201832_l0pt_aepcor.png");
  aa.SavePlot("el201832/m80to100/l0pt","norm xmax:100 suffix:_newtriggerSF_s3:dy save:el201832_l0pt_aepcor_tagpt45.png");
}

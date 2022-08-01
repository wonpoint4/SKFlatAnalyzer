//efficiency status report
#include "EfficiencyPlotter.cc"
void dirap(){
  EfficiencyPlotter aa;
  for(TString era:{"2016a","2016b","2017","2018"}){
    for(TString channel:{"ee","el"}){
      aa.SavePlot(channel+era+"/m80to100/leta","xmin:-2.5 xmax:2.5 norm rebin:2 xtite:'#eta(e)' preliminary save:eff_v7/"+channel+era+"_leta.png");
      aa.SavePlot(channel+era+"/m80to100/leta","xmin:-2.5 xmax:2.5 norm rebin:2 xtite:'#eta(e)' preliminary suffix:_noz0weight:sim save:eff_v7/"+channel+era+"_leta_noz0weight.png");
      aa.SavePlot(channel+era+"/m80to100/dirap","xmin:-2.4 xmax:2.4 norm rebin:2 xtite:'y(ee)' preliminary save:eff_v7/"+channel+era+"_dirap.png");
      aa.SavePlot(channel+era+"/m80to100/dirap","xmin:-2.4 xmax:2.4 norm rebin:2 xtite:'y(ee)' preliminary suffix:_noz0weight:sim save:eff_v7/"+channel+era+"_dirap_noz0weight.png");
      aa.SavePlot("v7_3/"+channel+era+"/m80to100/leta","xmin:-2.5 xmax:2.5 norm rebin:2 xtite:'#eta(e)' preliminary save:eff_v7_3/"+channel+era+"_leta.png");
      aa.SavePlot("v7_3/"+channel+era+"/m80to100/leta","xmin:-2.5 xmax:2.5 norm rebin:2 xtite:'#eta(e)' preliminary suffix:_noz0weight:sim save:eff_v7_3/"+channel+era+"_leta_noz0weight.png");
      aa.SavePlot("v7_3/"+channel+era+"/m80to100/dirap","xmin:-2.4 xmax:2.4 norm rebin:2 xtite:'y(ee)' preliminary save:eff_v7_3/"+channel+era+"_dirap.png");
      aa.SavePlot("v7_3/"+channel+era+"/m80to100/dirap","xmin:-2.4 xmax:2.4 norm rebin:2 xtite:'y(ee)' preliminary suffix:_noz0weight:sim save:eff_v7_3/"+channel+era+"_dirap_noz0weight.png");
    }
  }
}
void lpt(){
  EfficiencyPlotter aa;
  for(TString era:{"2016a","2016b","2017","2018"}){
    for(TString channel:{"ee","el"}){
      aa.SavePlot(channel+era+"/m80to100/lpt","xmax:100 norm xtite:'p_{T}(e) [GeV]' preliminary save:"+channel+era+"_lpt.png");
      aa.SavePlot(channel+era+"/m80to100/lpt","xmax:100 norm xtite:'p_{T}(e) [GeV]' preliminary suffix:_electronIDSF_s7_m0:sim save:"+channel+era+"_lpt_nofix.png");
      aa.SavePlot(channel+era+"/m80to100/lpt","xmax:100 norm xtite:'p_{T}(e) [GeV]' preliminary suffix:_noprefireweight:sim save:"+channel+era+"_lpt_noprefireweight.png");
    }
  }
}  
void plotall(){
  dirap();
  //lpt();
}

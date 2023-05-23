#include"EfficiencyPlotter.cc"
void plot(){
  EfficiencyPlotter aa;
  aa.plotdir="fig/220728";
  for(TString era:{"2016a","2016b"}){
    aa.SavePlot("ee"+era+"/m80to100/leta","norm xmin:-2.5 xmax:2.5 rebin:2 xtitle:'#eta(e)' save:ee"+era+"_leta");
    aa.SavePlot("fine/ee"+era+"/m80to100/leta","norm xmin:-2.5 xmax:2.5 rebin:2 xtitle:'#eta(e)' save:ee"+era+"_leta_fine");
    aa.SavePlot("ee"+era+"/m80to100/dirap","norm xmin:-2.4 xmax:2.4 rebin:2 xtitle:'y(ee)' save:ee"+era+"_rapidity");
    aa.SavePlot("fine/ee"+era+"/m80to100/dirap","norm xmin:-2.4 xmax:2.4 rebin:2 xtitle:'y(ee)' save:ee"+era+"_rapidity_fine");
    aa.SavePlot("el"+era+"/m80to100/leta","norm xmin:-2.5 xmax:2.5 rebin:2 xtitle:'#eta(e)' save:el"+era+"_leta");
    aa.SavePlot("fine/el"+era+"/m80to100/leta","norm xmin:-2.5 xmax:2.5 rebin:2 xtitle:'#eta(e)' save:el"+era+"_leta_fine");
    aa.SavePlot("el"+era+"/m80to100/dirap","norm xmin:-2.4 xmax:2.4 rebin:2 xtitle:'y(ee)' save:el"+era+"_rapidity");
    aa.SavePlot("fine/el"+era+"/m80to100/dirap","norm xmin:-2.4 xmax:2.4 rebin:2 xtitle:'y(ee)' save:el"+era+"_rapidity_fine");
  }
}

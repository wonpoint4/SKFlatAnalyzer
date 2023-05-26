#include"EfficiencyPlotter.cc"
void plot(){
  EfficiencyPlotter aa;
  aa.plotdir="fig/230524";
  for(TString channel:{"ee","el","mm","mu"}){
    for(TString era:{"2016a","2016b","2017","2018"}){
      TString common=" norm chi2detail pvalue shapesys 2:finey sysname:effAN2 ";
      if(channel.BeginsWith('e')){
	aa.SavePlotCondor(channel+era+"/m80to100/l1eta",common+" xmin:-2.5 xmax:2.5 rebin:2 xtitle:'electron #eta' save:"+channel+era+"_l1eta.png ");
	aa.SavePlotCondor(channel+era+"/m80to100/l1pt",common+" xmin:15 xmax:100 xtitle:'electron p_{T} [GeV]' rebin:{10,15,20,25,30,35,40,45,50,60,70,100} save:"+channel+era+"_l1pt.png ");
      }else{
	aa.SavePlotCondor(channel+era+"/m80to100/l1eta",common+" xmin:-2.4 xmax:2.4 rebin:2 xtitle:'muon #eta' save:"+channel+era+"_l1eta.png ");
	aa.SavePlotCondor(channel+era+"/m80to100/l1pt",common+" xmin:10 xmax:100 xtitle:'muon p_{T} [GeV]' rebin:{10,15,20,25,30,35,40,45,50,70,100} save:"+channel+era+"_l1pt.png ");
      }
      aa.SavePlotCondor(channel+era+"/m80to100/dimass",common+" xtitle:'m(ll) [GeV]' save:"+channel+era+"_dimass.png ");
      aa.SavePlotCondor(channel+era+"/m80to100/dirap",common+" xmin:-2.4 xmax:2.4 rebin:2 xtitle:'y(ll)' save:"+channel+era+"_dirap.png ");
      aa.SavePlotCondor(channel+era+"/m80to100/dipt",common+" xmax:100 rebin:2 xtitle:'p_{T}(ll) [GeV]' save:"+channel+era+"_dipt.png ");
      aa.SavePlotCondor(channel+era+"/m80to100/cost",common+" absx xtitle:'cos(#theta_{CS})' save:"+channel+era+"_cost.png ");
    }
  }
}

#include"EfficiencyPlotter.cc"

void plot(){
  EfficiencyPlotter aa;
  aa.plotdir="fig/230320eff_v13_v12/";
  for(TString channel:{"ee","el"}){
    for(TString era:{"2016a","2016b","2017","2018","201727","201732","201828","201832"}){
      if(channel=="ee"){
	if(era.Length()==6) continue;
      }
      aa.SavePlotCondor(channel+era+"/m80to100/l0pt","xmax:100 norm sysname:efficiencySF save:"+channel+era+"_l0pt_v13");
      aa.SavePlotCondor(channel+era+"/m80to100/l1pt","xmax:100 norm sysname:efficiencySF save:"+channel+era+"_l1pt_v13");
      aa.SavePlotCondor(channel+era+"/m80to100/lpt","xmax:100 norm sysname:efficiencySF save:"+channel+era+"_lpt_v13");
      aa.SavePlotCondor("v12/"+channel+era+"/m80to100/l0pt","xmax:100 norm sysname:efficiencySF save:"+channel+era+"_l0pt_v12");
      aa.SavePlotCondor("v12/"+channel+era+"/m80to100/l1pt","xmax:100 norm sysname:efficiencySF save:"+channel+era+"_l1pt_v12");
      aa.SavePlotCondor("v12/"+channel+era+"/m80to100/lpt","xmax:100 norm sysname:efficiencySF save:"+channel+era+"_lpt_v12");
    }
  }
}

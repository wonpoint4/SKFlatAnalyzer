#include"EfficiencyPlotter.cc"

void plot1(){
  EfficiencyPlotter aa;
  //aa.entries[0].styles[0].markersize=0.4;
  aa.plotdir="fig/230625_2/costym";
  //TString common=" xmin:0 xmax:660 chi2detail pvalue norm decotest ";
  TString common=" chi2detail pvalue norm decotest Zmin:82 Zmax:100 decoeffstat ";
  for(TString era:{"2016a","2016b","2017","2018","Run2"}){
    TString era_= era=="Run2"? "201[678][ab]?" : era;
    for(TString channel:{"ee","el","mm","mu"}){
      aa.SavePlotCondor(channel+era_+"/costym","sysname:effAN2 save:"+channel+era+common);
      aa.SavePlotCondor("v17/"+channel+era_+"/costym","sysname:effAN2new save:"+channel+era+"_v17"+common);
      if(channel[0]=='e')
	aa.SavePlotCondor("v17/"+channel+era_+"/costym","sysname:effAN2newpp save:"+channel+era+"_v17_pp"+common);
    }
  }
}
void plot2(){
  EfficiencyPlotter aa;
  aa.plotdir="fig/230625_2/ym";
  TString common=" chi2detail pvalue norm decotest project:yz ";
  for(TString era:{"2016a","2016b","2017","2018","Run2"}){
    TString era_= era=="Run2"? "201[678][ab]?" : era;
    for(TString channel:{"ee","el","mm","mu"}){
      aa.SavePlotCondor(channel+era_+"/costym","sysname:effAN2 save:"+channel+era+"_deco_noresidual"+common);
      if(channel[0]=='e')
	aa.SavePlotCondor(channel+era_+"/costym","suffix:_electronIDSF_s12_m0:dy save:"+channel+era+"_deco_residual"+common);
      if(channel[0]=='m')
	aa.SavePlotCondor(channel+era_+"/costym","suffix:_muonIDSF_s11_m0:dy save:"+channel+era+"_deco_residual"+common);
    }
  }
}
void plot3(){
  EfficiencyPlotter aa;
  aa.plotdir="fig/230625_2/costy";
  TString common=" Zmin:82 Zmax:100 chi2detail pvalue norm project:xy decotest decoeffstat ";
  for(TString era:{"2016a","2016b","2017","2018","Run2"}){
    TString era_= era=="Run2"? "201[678][ab]?" : era;
    for(TString channel:{"ee","el","mm","mu"}){
      aa.SavePlotCondor(channel+era_+"/costym"," sysname:effAN2 save:"+channel+era+"_costy"+common);
      aa.SavePlotCondor("v17/"+channel+era_+"/costym"," sysname:effAN2new save:"+channel+era+"_costy_v17"+common);
      if(channel[0]=='e')
	aa.SavePlotCondor("v17/"+channel+era_+"/costym"," sysname:effAN2newpp save:"+channel+era+"_costy_v17_pp"+common);
    }
  }
}
void plot4(){
  EfficiencyPlotter aa;
  //aa.entries[0].styles[0].markersize=0.4;
  aa.plotdir="fig/230625_2/letapt";
  for(TString era:{"2016a","2016b","2017","2018","Run2"}){
    TString era_= era=="Run2"? "201[678][ab]?" : era;
    TString common=" chi2 pvalue norm decoeffstat ";
    for(TString channel:{"ee","el","mm","mu"}){
      aa.SavePlotCondor(channel+era_+"/m80to100/letapt","sysname:effAN2 chi2detail save:"+channel+era+"_letapt"+common);
      aa.SavePlotCondor("v17/"+channel+era_+"/m80to100/letapt","sysname:effAN2new chi2detail save:"+channel+era+"_letapt_v17"+common);
      if(channel[0]=='e')
	aa.SavePlotCondor("v17/"+channel+era_+"/m80to100/letapt","sysname:effAN2newpp chi2detail save:"+channel+era+"_letapt_v17_pp"+common);
    }
  }
}
void plot1d(){
  EfficiencyPlotter aa;
  //aa.entries[0].styles[0].markersize=0.4;
  aa.plotdir="fig/230625_2/1d";
  for(TString era:{"2016a","2016b","2017","2018","Run2"}){
    TString era_= era=="Run2"? "201[678][ab]?" : era;
    for(TString channel:{"ee","el","mm","mu"}){
      TString common=" chi2 pvalue norm decoeffstat ";
      if(channel[0]=='e'){
	aa.SavePlotCondor(channel+era_+"/m80to100/letapt","sysname:effAN2 chi2detail project:x xmin:-2.5 xmax:2.5 save:"+channel+era+"_leta"+common);
	aa.SavePlotCondor(channel+era_+"/m80to100/letapt","sysname:effAN2 chi2detail project:y xmin:15 xmax:500 save:"+channel+era+"_lpt"+common);
	aa.SavePlotCondor("v17/"+channel+era_+"/m80to100/letapt","sysname:effAN2new chi2detail project:x xmin:-2.5 xmax:2.5 save:"+channel+era+"_leta_v17"+common);
	aa.SavePlotCondor("v17/"+channel+era_+"/m80to100/letapt","sysname:effAN2new chi2detail project:y xmin:15 xmax:500 save:"+channel+era+"_lpt_v17"+common);
	aa.SavePlotCondor("v17/"+channel+era_+"/m80to100/letapt","sysname:effAN2newpp chi2detail project:x xmin:-2.5 xmax:2.5 save:"+channel+era+"_leta_v17_pp"+common);
	aa.SavePlotCondor("v17/"+channel+era_+"/m80to100/letapt","sysname:effAN2newpp chi2detail project:y xmin:15 xmax:500 save:"+channel+era+"_lpt_v17_pp"+common);
      }else{
	aa.SavePlotCondor(channel+era_+"/m80to100/letapt","sysname:effAN2 chi2detail project:x xmin:-2.4 xmax:2.4 save:"+channel+era+"_leta"+common);
	aa.SavePlotCondor(channel+era_+"/m80to100/letapt","sysname:effAN2 chi2detail project:y xmin:10 xmax:200 save:"+channel+era+"_lpt"+common);
	aa.SavePlotCondor("v17/"+channel+era_+"/m80to100/letapt","sysname:effAN2new chi2detail project:x xmin:-2.4 xmax:2.4 save:"+channel+era+"_leta_v17"+common);
	aa.SavePlotCondor("v17/"+channel+era_+"/m80to100/letapt","sysname:effAN2new chi2detail project:y xmin:10 xmax:200 save:"+channel+era+"_lpt_v17"+common);
      }	
      common=" Zmin:82 Zmax:100 chi2detail pvalue norm ";
      aa.SavePlotCondor(channel+era_+"/costym"," sysname:effAN2 project:x save:"+channel+era+"_cost"+common);
      aa.SavePlotCondor("v17/"+channel+era_+"/costym"," sysname:effAN2new project:x save:"+channel+era+"_cost_v17"+common);
      if(channel[0]=='e'){
	aa.SavePlotCondor("v17/"+channel+era_+"/costym"," sysname:effAN2newpp project:x save:"+channel+era+"_cost_v17_pp"+common);
	aa.SavePlotCondor("v17/"+channel+era_+"/costym"," sysname:effAN2newpp project:y xmin:-2.4 xmax:2.4 save:"+channel+era+"_y_v17_pp"+common);
      }
      common=" chi2detail pvalue norm rebin:2";
      aa.SavePlotCondor(channel+era_+"/m80to100/dirap"," sysname:effAN2 xmin:-2.4 xmax:2.4 save:"+channel+era+"_y"+common);
      aa.SavePlotCondor("v17/"+channel+era_+"/m80to100/dirap"," sysname:effAN2new xmin:-2.4 xmax:2.4 save:"+channel+era+"_y_v17"+common);
      if(channel[0]=='e'){
	aa.SavePlotCondor("v17/"+channel+era_+"/m80to100/dirap"," sysname:effAN2newpp xmin:-2.4 xmax:2.4 save:"+channel+era+"_y_v17_pp"+common);
      }
    }
  }  
}
void plot(){
  plot1();
  plot3();
  plot4();
  plot1d();
}

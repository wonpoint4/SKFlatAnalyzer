#include"EfficiencyPlotter.cc"

void plot1(){
  EfficiencyPlotter aa;
  //aa.entries[0].styles[0].markersize=0.4;
  aa.plotdir="fig/230620/costym_test";
  //TString common=" xmin:0 xmax:660 chi2detail pvalue norm decotest ";
  TString common=" chi2detail pvalue norm decotest Zmin:82 Zmax:100 ";
  for(TString era:{"2016a","2016b","2017","2018","Run2"}){
    TString era_= era=="Run2"? "201[678][ab]?" : era;
    for(TString channel:{"ee","el","mm","mu"}){
      aa.SavePlotCondor(channel+era_+"/costym","sysname:effAN3 save:"+channel+era+"_deco"+common);
      aa.SavePlotCondor(channel+era_+"/costym","sysname:effAN2 save:"+channel+era+"_deco_noresidual"+common);
      if(channel[0]=='e')
	aa.SavePlotCondor(channel+era_+"/costym","suffix:_electronIDSF_s12_m0:dy save:"+channel+era+"_deco_residual"+common);
      if(channel[0]=='m')
	aa.SavePlotCondor(channel+era_+"/costym","suffix:_muonIDSF_s11_m0:dy save:"+channel+era+"_deco_residual"+common);
    }
  }
}
void plot2(){
  EfficiencyPlotter aa;
  aa.plotdir="fig/230620/ym";
  TString common=" chi2detail pvalue norm decotest project:yz ";
  for(TString era:{"2016a","2016b","2017","2018","Run2"}){
    TString era_= era=="Run2"? "201[678][ab]?" : era;
    for(TString channel:{"ee","el","mm","mu"}){
      aa.SavePlotCondor(channel+era_+"/costym","sysname:effAN3 save:"+channel+era+"_deco"+common);
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
  aa.plotdir="fig/230620/costy";
  TString common=" Zmin:82 Zmax:100 chi2detail pvalue norm project:xy decotest ";
  for(TString era:{"2016a","2016b","2017","2018","Run2"}){
    TString era_= era=="Run2"? "201[678][ab]?" : era;
    for(TString channel:{"ee","el","mm","mu"}){
      aa.SavePlotCondor(channel+era_+"/costym"," sysname:effAN2 save:"+channel+era+"_costy"+common);
      if(channel[0]=='e'){
	aa.SavePlotCondor("v16/"+channel+era_+"/costym"," sysname:effAN2new save:"+channel+era+"_costy_v16"+common);
	aa.SavePlotCondor("v16/"+channel+era_+"/costym"," sysname:effAN2newpp save:"+channel+era+"_costy_v16_pp"+common);
      }
    }
  }
}
void plot4(){
  EfficiencyPlotter aa;
  //aa.entries[0].styles[0].markersize=0.4;
  aa.plotdir="fig/230620/letapt";
  for(TString era:{"2016a","2016b","2017","2018","Run2"}){
    TString era_= era=="Run2"? "201[678][ab]?" : era;
    TString common=" chi2 pvalue norm decoeffstat";
    for(TString channel:{"ee","el","mm","mu"}){
      aa.SavePlotCondor(channel+era_+"/m80to100/letapt","sysname:effAN2 chi2detail save:"+channel+era+"_letapt"+common);
      if(channel[0]=='e'){
	aa.SavePlotCondor("v16/"+channel+era_+"/m80to100/letapt","sysname:effAN2new chi2detail save:"+channel+era+"_letapt_v16"+common);
	aa.SavePlotCondor("v16/"+channel+era_+"/m80to100/letapt","sysname:effAN2newpp chi2detail save:"+channel+era+"_letapt_v16_pp"+common);
      }
    }
  }
}
void plot(){
  plot3();
  plot4();
}

#include"Chi2Plotter.cc"

void save_compare(TString suffix=""){
  Chi2Plotter aa;
  TString sysname="totalsys"+suffix;
  aa.plotdir="fig/240125/";
  for(TString channel:{"el","mu","ee","mm"}){
    for(TString era:{"2016a","2016b","2017","2018","Run2"}){
      TString era_=era;
      if(era=="Run2") era_="201[678][ab]?";
      TString common=" chi2 movesys ";

      aa.SavePlotCondor(channel+era_+"/ym","sysname:"+sysname+" save:"+channel+era+"_ym"+suffix+".png "+common);
      aa.SavePlotCondor(channel+era_+"/ym","sysname:"+sysname+" project:x save:"+channel+era+"_y"+suffix+".png "+common);
      aa.SavePlotCondor(channel+era_+"/ym","sysname:"+sysname+" widthweight project:y save:"+channel+era+"_m"+suffix+".png "+common);
      aa.SavePlotCondor(channel+era_+"/l0etapt","sysname:"+sysname+" save:"+channel+era+"_l0etapt"+suffix+".png "+common);
      aa.SavePlotCondor(channel+era_+"/l1etapt","sysname:"+sysname+" save:"+channel+era+"_l1etapt"+suffix+".png "+common);
      aa.SavePlotCondor(channel+era_+"/l0etapt","sysname:"+sysname+" project:x save:"+channel+era+"_l0eta"+suffix+".png "+common);
      aa.SavePlotCondor(channel+era_+"/l0etapt","sysname:"+sysname+" logx widthweight project:y save:"+channel+era+"_l0pt"+suffix+".png "+common);
    }
  }
}

void save(){
  save_compare("");
  save_compare("_noresidual");
  save_compare("_noresidual_correlated");
  save_compare("_correlated");
}

#include"AFBSystPlotter.cc"

void save_compare(TString suffix="_more"){
  AFBSystPlotter aa;
  TString sysname="totalsys"+suffix;
  aa.plotdir="fig/231217/";
  for(TString channel:{"ee","mm"}){
    for(TString era:{"2016a","2016b","2017","2018","No2016a","Run2"}){
      for(TString region:{"0bjet","nbjet"}){
	TString era_=era;
	if(era=="Run2") era_="201[678][ab]?";
	if(era=="No2016a") era_="201[678][b]?";

	TString common=" chi2 movesys";
      
	aa.SavePlotCondor(channel+era_+"/"+region+"/dimass","widthweight logx xmin:52 xmax:3000 sysname:"+sysname+" save:"+channel+era+"_"+region+"_dimass"+suffix+".png "+common);
	aa.SavePlotCondor(channel+era_+"/"+region+"/dirap","widthweight xmin:-2.4 xmax:2.4 sysname:"+sysname+" save:"+channel+era+"_"+region+"_dirap"+suffix+".png "+common);
	aa.SavePlotCondor(channel+era_+"/"+region+"/dipt","widthweight logx xmin:2 xmax:650 sysname:"+sysname+" save:"+channel+era+"_"+region+"_dipt"+suffix+".png "+common);
	for(int im=0;im<4;im++){
	  aa.SavePlotCondor(channel+era_+"/"+region+Form("/dirap_m%d",im),"widthweight xmin:-2.4 xmax:2.4 sysname:"+sysname+" save:"+channel+era+"_"+region+Form("_dirap_m%d%s.png",im,suffix.Data())+common);
	  aa.SavePlotCondor(channel+era_+"/"+region+Form("/dipt_m%d",im),"widthweight logx xmin:2 xmax:650 sysname:"+sysname+" save:"+channel+era+"_"+region+Form("_dipt_m%d%s.png",im,suffix.Data())+common);
	}
      }
    }
  }
}
void save_diff(TString suffix="_more"){
  AFBSystPlotter aa("data-mi-tau_mi-vv-wjets-tt-st-qcdss-aa");
  TString sysname="totalsys"+suffix;
  aa.plotdir="fig/231217/";
  for(TString channel:{"ee","mm"}){
    for(TString era:{"2016a","2016b","2017","2018","No2016a","Run2"}){
      for(TString region:{"0bjet","nbjet"}){
	TString era_=era;
	if(era=="Run2") era_="201[678][ab]?";
	if(era=="No2016a") era_="201[678][b]?";
      
	TString common=" chi2";
	aa.SavePlotCondor(channel+era_+"/"+region+"/dimass","widthweight logx xmin:52 xmax:3000 sysname:"+sysname+" save:"+channel+era+"_"+region+"_dimass"+suffix+"_diff.png "+common);
	aa.SavePlotCondor(channel+era_+"/"+region+"/dirap","widthweight xmin:-2.4 xmax:2.4 sysname:"+sysname+" save:"+channel+era+"_"+region+"_dirap"+suffix+"_diff.png "+common);
	aa.SavePlotCondor(channel+era_+"/"+region+"/dipt","widthweight logx xmin:2 xmax:650 sysname:"+sysname+" save:"+channel+era+"_"+region+"_dipt"+suffix+"_diff.png "+common);
	for(int im=0;im<4;im++){
	  aa.SavePlotCondor(channel+era_+"/"+region+Form("/dirap_m%d",im),"widthweight xmin:-2.4 xmax:2.4 sysname:"+sysname+" save:"+channel+era+"_"+region+Form("_dirap_m%d%s_diff.png",im,suffix.Data())+common);
	  aa.SavePlotCondor(channel+era_+"/"+region+Form("/dipt_m%d",im),"widthweight logx xmin:2 xmax:650 sysname:"+sysname+" save:"+channel+era+"_"+region+Form("_dipt_m%d%s_diff.png",im,suffix.Data())+common);
	}
      }
    }
  }
}

void save(){
  save_compare("_more2");
}

#include"AFBPlotter.cc"
void plot(){
  AFBPlotter aa;
  aa.plotdir="fig/230725";
  for(TString channel:{"ee","mm"}){
    for(TString era:{"2016a","2016b","2017","2018","Run2"}){
      for(TString region:{"0bjet","nbjet"}){
	TString era_= era=="Run2" ? "201[678][ab]?" : era;
	aa.SavePlotCondor(channel+era_+"/"+region+"/jets","Xmin:52 Xmax:3000 rebinU:{0,10} absY chi2detail pvalue sysname:totalsys norm xdivision:608 save:"+channel+era+"_"+region);
      }
    }
  }
}
